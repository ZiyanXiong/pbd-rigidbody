classdef Model < handle
	%Model Class to hold model data

	properties
		bodies % list of bodies
		constraints %list of constraints
		collider % collision handler
        constraintList % lise of constraints and contacts 
		grav % gravity
		ground % ground transform, with Z up
        bodyLayers % Bodyies at each layer
        constraintLayers% constraints indices at each layer
		
		t % current time
		h % time step
		tEnd % end time
		k % current step

		steps % number of steps to take (tEnd/h)
		substeps % number of substeps per step
		iters % number of iterations per substep
		ks % current substep
		hs % sub time step

		% Energy
		tt % time
		TT % kinetic energy
		VV % potential energy
		plotH % whether to plot the energy at the end
		computeH % whether to compute the energy
		Hexpected % expected energy

		% Drawing etc. (not required for sim)
		name % scene name
		drawHz % refresh rate (0 for no draw)
		view % initial viewing angle
		axis % initial axis
		video %
	end

	methods
		function this = Model()
			apbd.Model.clearGlobal();
			
			% Default values
			this.bodies = {};
			this.constraints = {};
			this.collider = apbd.Collider(this);
			this.grav = [0 0 -980]';
			this.ground.E = zeros(4);
            this.constraintLayers = {};

			this.t = 0;
			this.h = 1/30;
			this.tEnd = 1;
			this.k = 0;

			this.steps = 0; % will be computed in init()
			this.substeps = 10;
			this.iters = 1;
			this.ks = 0;
			this.hs = this.h/this.substeps;

			this.tt = [];
			this.TT = [];
			this.VV = [];
			this.plotH = false;
			this.computeH = true;
			this.Hexpected = zeros(1,2);

			% Drawing etc.
			this.name = 'Default';
			this.ground.size = 10;
			this.drawHz = 15;
			this.view = 3;
			this.axis = [];
			this.video = [];
		end

		%%
		function init(this)
			if usejava('jvm')
				colormap('default'); % Restore colormap mode
			end

			for i = 1 : length(this.bodies)
				this.bodies{i}.init();
			end
			for i = 1 : length(this.constraints)
				this.constraints{i}.init();
			end

			% Other initial values
			this.steps = ceil(this.tEnd/this.h);
			this.hs = this.h/this.substeps;
			this.k = 0;
			this.ks = 0;
			this.computeEnergies();
			this.draw();
		end

		%%
		function simulate(this)
			while this.k < this.steps
				this.ks = 0;
                this.clearBodyShockPropInfo();
                this.collider.run();
                this.constructConstraintGraph();
				while this.ks < this.substeps
                    this.stepBDF1();
				    %this.solveConJacobi();
                    this.solveConSP();
                    %this.solveConGSSP();
                    this.solveConGS();
					this.t = this.t + this.hs;
					this.ks = this.ks + 1;
                    %this.rigidify();
				end
				this.k = this.k + 1;
				this.computeEnergies();
				this.draw();
			end
        end

        %%
        function applyimplulse(this, depth)
			this.draw();
            
            for i = 1:length(this.bodies)
                this.bodies{i}.x1 = this.bodies{i}.x;
                if(i ~= 1)
                    %this.bodies{i}.x = this.bodies{i -1}.x;
                end
                this.collider.collisions{end+1} = apbd.ConCollGroundRigid2dLinear(this.bodies{i},this.ground.E);
			    this.collider.collisions{end}.d = depth * i;
		        this.collider.collisions{end}.xl = [0.5 -0.5 0]';
			    this.collider.collisions{end}.nw = [0 1 0]';
                this.collider.collisions{end}.solveNorPos();
                this.bodies{i}.applyJacobi();
			    this.draw();
            end

            %{
			this.collider.collisions{end}.d = -depth;
            this.collider.collisions{end}.solveNorPos();
	        for i = 1 : length(this.bodies)
                this.bodies{i}.applyJacobi();
            end
			this.draw();
            %}
        end

        %%
        function c = computeC(this, p1, p2)
            nw = [0 1 0]';
            rl1 = [-0.5 -0.5 0]';
            rl2 = [0.0 -0.5 0]';
		    m1 = this.bodies{1}.Mp;
		    I1 = this.bodies{1}.Mr;

	        % Position update
	        dpw = p1*nw;
	        dp1 = dpw/m1;
	        % Quaternion update
            q1 = apbd.BodyRigid2d.unproj(this.bodies{1}.x1_0);
	        dpl1 = se3.qRotInv(q1,dpw);
	        qtmp1 = [se3.qRot(q1, I1.\se3.cross(rl1,dpl1)); 0];
            %qtmp1 = [I1.\se3.cross(rl1,dpl1); 0];
	        %dq1 = se3.qMul(sin(0.5 * qtmp1),q1);
            dq1 = 0.5 * se3.qMul(qtmp1,q1);
            dx = [0 0 0 0]';
            dx(1:2) = dq1(3:4);
            dx(3:4) = dp1(1:2);
            this.bodies{1}.x1 = this.bodies{1}.x1_0 + dx;
            
	        dpw = p2*nw;
	        dp1 = dpw/m1;
	        % Quaternion update
            q1 = apbd.BodyRigid2d.unproj(this.bodies{1}.x1_0);
	        dpl1 = se3.qRotInv(q1,dpw);
	        qtmp1 = [se3.qRot(q1, I1.\se3.cross(rl2,dpl1)); 0];
            %qtmp1 = [I1.\se3.cross(rl1,dpl1); 0];
	        %dq1 = se3.qMul(sin(0.5 * qtmp1),q1);
            dq1 = 0.5 * se3.qMul(qtmp1,q1);
            dx = [0 0 0 0]';
            dx(1:2) = dq1(3:4);
            dx(3:4) = dp1(1:2);
            this.bodies{1}.x1 = this.bodies{1}.x1 + dx;
            
            this.rigidify();

	        xw1 = this.bodies{1}.transformPoint(rl1);
	        [q, p] = apbd.BodyRigid2d.unproj(this.bodies{1}.x0);
	        xwi = se3.qRot(q,rl1) + p;
            d1 = dot(nw, xw1 - xwi);
            xw2 = this.bodies{1}.transformPoint(rl2);
            xwi = se3.qRot(q,rl2) + p;
            d2 = dot(nw, xw2 - xwi);
            c = 0.5 * d1^2 + 0.5 * d2^2;
        end

        %%
        function plotC(this)
            %this.bodies{1}.x = this.bodies{1}.x +  [0.2 0 0 -0.02]';
            this.bodies{1}.x1 = this.bodies{1}.x;
            this.bodies{1}.x1_0 = this.bodies{1}.x1;
            this.bodies{1}.x0 = this.bodies{1}.x + [0 0 0 0.5]';
            
            this.collider.collisions{end+1} = apbd.ConCollGroundRigid2d(this.bodies{1},this.ground.E);
		    this.collider.collisions{end}.d = 0;
	        this.collider.collisions{end}.xl = [-0.5 -0.5 0]';
		    this.collider.collisions{end}.nw = [0 1 0]';

            this.collider.collisions{end+1} = apbd.ConCollGroundRigid2d(this.bodies{1},this.ground.E);
		    this.collider.collisions{end}.d = 0;
	        this.collider.collisions{end}.xl = [0.0 -0.5 0]';
		    this.collider.collisions{end}.nw = [0 1 0]';

            % Plot the value of C
            p1 = -0.1:0.02:1.2;
            p2 = -0.1:0.02:1.2;
            [p1g, p2g] = meshgrid(p1,p2) ;
            c = zeros(length(p1));
            
            for i = 1:length(p1)
                for j = 1:length(p2)
                    c(i,j) = this.computeC(p1g(i,j), p2g(i,j));
                end
            end

            surf(p1g,p2g,c)
            hold on;

            p1_history = zeros(this.iters, 1);
            p2_history = zeros(this.iters, 1);
            c_history = zeros(this.iters, 1);
            c_acc = zeros(this.iters, 1);

            this.bodies{1}.x = this.bodies{1}.x1_0;
            this.bodies{1}.x1 = this.bodies{1}.x1_0;
			for i = 1 : length(this.constraints)
				this.constraints{i}.clear();
            end
		    %this.collider.run();
            collisions = this.collider.collisions;
			for iter = 1 : this.iters % index 0
				%fprintf('  iter %d\n',iter);
				% Clear the Jacobi updates
				for i = 1 : length(this.bodies)
					this.bodies{i}.clearJacobi();
                end

				for j = 1 : length(collisions)
					collisions{j}.update();
                end

                p1_history(iter) = collisions{1}.lambda(1);
                p2_history(iter) = collisions{2}.lambda(1);
                c_history(iter) = 0.5 * collisions{1}.d^2 + 0.5 * collisions{2}.d^2;

				% Gauss-Seidel solve for non-collision constraints
				for j = 1 : length(this.constraints)
					this.constraints{j}.solve();
                end
				% Update the collisions with the latest body states
				for j = 1 : length(collisions)
					collisions{j}.update();
                    collisions{j}.solveNorPos();
                    %collisions{j}.solveTanVel(this.k,this.ks,this.hs);
    		        for i = 1 : length(this.bodies)
                        this.bodies{i}.applyJacobi();
                    end
                end

		        for i = 1 : length(this.bodies)
                    %this.bodies{i}.applyJacobi();
                end
            end

            this.bodies{1}.x = this.bodies{1}.x1_0;
            this.bodies{1}.x1 = this.bodies{1}.x1_0;
            for iter = 1 : this.iters % index 0
                c_acc(iter) = this.computeC(p1_history(iter), p2_history(iter));
            end
            plot3(p1_history, p2_history, c_history,'-*', color="red")
            hold on;
            %plot3(p1_history, p2_history, c_acc,'-*', color="yellow")
        end

		%%
		function stepBDF1(this)
			for i = 1 : length(this.bodies)
				this.bodies{i}.stepBDF1(this.k,this.ks,this.hs,this.grav);
			end
        end

		%%
        function rigidify(this)
			for i = 1 : length(this.bodies)
                this.bodies{i}.regularize();
			end
        end

		%%
        function clearBodyShockPropInfo(this)
			for i = 1 : length(this.bodies)
                this.bodies{i}.shockParentConIndex = {};
                this.bodies{i}.conIndex = [];
                this.bodies{i}.layer = 99;
            end
            this.bodyLayers = {};
            this.constraintLayers = {};
        end

		%%
        function constructConstraintGraph(this)
            this.constraintList = cat(2, this.constraints, this.collider.collisions);
            bodylayer = [];
            conlayer = [];
			for i = 1 : length(this.constraintList)
                if this.constraintList{i}.ground == true
                    bodylayer(end+1) = this.constraintList{i}.bodies{1}.index;
                    this.constraintList{i}.bodies{1}.layer = 1;
                    this.constraintList{i}.bodies{1}.shockParentConIndex{end + 1} = i;
                    this.constraintList{i}.bodies{1}.conIndex(end+1) = i;
                    conlayer(end+1) = i;
                end

                if length(this.constraintList{i}.bodies) == 2
                    this.constraintList{i}.bodies{1}.conIndex(end+1) = i;
                    this.constraintList{i}.bodies{2}.conIndex(end+1) = i;
                end
            end
            bodylayer = unique(bodylayer);
            conlayer = unique(conlayer);
            this.bodyLayers{end + 1} = bodylayer;
            this.constraintLayers{end + 1} = conlayer;

            for i = 1: length(this.bodies)
                this.bodies{i}.conIndex = unique(this.bodies{i}.conIndex);
            end
            
            layer = 2;
            while true
                bodylayer = [];
                conlayer = [];
                for i = 1 : length(this.bodyLayers{layer-1})
                     bodyIndex = this.bodyLayers{layer-1}(i);
                    for j = 1 : length(this.bodies{bodyIndex}.conIndex)
                        conIndex = this.bodies{bodyIndex}.conIndex(j);
                        if length(this.constraintList{conIndex}.bodies) == 1
                            continue;
                        end
                        if this.constraintList{conIndex}.bodies{1}.layer == this.constraintList{conIndex}.bodies{2}.layer
                            this.constraintLayers{layer-1}(end+1) = conIndex;
                            continue;
                        elseif this.constraintList{conIndex}.bodies{1}.layer > this.constraintList{conIndex}.bodies{2}.layer
                            %temp = this.constraintList{conIndex}.bodies{1};
                            %this.constraintList{conIndex}.bodies{1} = this.constraintList{conIndex}.bodies{2};
                            %this.constraintList{conIndex}.bodies{2} = temp;
                            if ismethod(this.constraintList{conIndex},'swapBody')
                                this.constraintList{conIndex}.swapBody();
                            end
                        end
                        if this.bodies{bodyIndex} == this.constraintList{conIndex}.bodies{2}
                            continue;
                        end
                        body2 = this.constraintList{conIndex}.bodies{2};
                        body2.shockParentConIndex{end + 1} = conIndex;
                        body2.layer = layer;
                        bodylayer(end+1) = body2.index;
                        conlayer(end+1) = conIndex;
                    end
                end
                if isempty(bodylayer)
                    break;
                end
                bodylayer = unique(bodylayer);
                conlayer = unique(conlayer);
                this.bodyLayers{end + 1} = bodylayer;
                this.constraintLayers{end + 1} = conlayer;
                layer = layer + 1;
            end

		end

		%%
		function solveConJacobi(this)
			%this.draw();
			%fprintf('substep %d\n',this.ks);
			for i = 1 : length(this.constraints)
				this.constraints{i}.clear();
            end
		    %this.collider.run();
            collisions = this.collider.collisions;
			for iter = 0 : this.iters-1 % index 0
				%fprintf('  iter %d\n',iter);
				% Clear the Jacobi updates
				for i = 1 : length(this.bodies)
					this.bodies{i}.clearJacobi();
				end
				% Gauss-Seidel solve for non-collision constraints
				for j = 1 : length(this.constraints)
					this.constraints{j}.solve();
				end

				% Solve all collision normals at the position level
				%fprintf('    ');
			    %this.collider.run();
                %collisions = this.collider.collisions;

				% Update the collisions with the latest body states
				for j = 1 : length(collisions)
					collisions{j}.update();
                end
				for j = 1 : length(collisions)
					    collisions{j}.solveNorPos();
				end
				% Apply Jacobi updates to the bodies
				%fprintf('    ');
				for i = 1 : length(this.bodies)
					this.bodies{i}.applyJacobi();
                end
                this.draw()

				% Update the collisions with the latest body states
				for j = 1 : length(collisions)
					collisions{j}.update();
                end
                
			    % Solve all collision tangents at the velocity level
			    for j = 1 : length(collisions)
				    collisions{j}.solveTanVel(this.k,this.ks,this.hs);
				    %fprintf('%e\n',collisions{j}.C(1));
			    end
			    % Apply Jacobi updates
			    %fprintf('    ');
			    for i = 1 : length(this.bodies)
				    this.bodies{i}.applyJacobi();
				    %fprintf('%e ',this.bodies{i}.x(end));
			    end
			    %fprintf('\n');
                %this.draw()
            end
		end

		%%
        function solveConSP(this)
			%this.draw();
			%fprintf('substep %d\n',this.ks);
			for i = 1 : length(this.constraints)
				this.constraints{i}.clear();
            end
            
			for i = 1 : length(this.bodies)
                for j = 1 : length(this.bodies{i}.shockParentConIndex)
				    conIndex = this.bodies{i}.shockParentConIndex{j};
                    this.constraintList{conIndex}.shockProp = true;
                end
            end
            
            for i = 1 : length(this.constraintLayers)
                for iter = 1:this.iters
                    for j = 1 : length(this.constraintLayers{i})
                        this.constraintList{this.constraintLayers{i}(j)}.solve(this.hs);
                    end
                end
            end
            
            %{
            for i = 1:length(this.bodyLayers)
                for iter = 1:50
                    for j = 1:length(this.bodyLayers{i})
                        body = this.bodies{this.bodyLayers{i}(j)};
                        for k = 1:length(body.shockParentConIndex)
                            %conIndex = body.shockParentConIndex{k};
                            %this.constraintList{conIndex}.shockProp = true;
                        end
                        for k = 1:length(body.conIndex)
                            conIndex = body.conIndex(k);
                            this.constraintList{conIndex}.solve(this.hs);
                        end
                        
                    end
                end
            end 
            %}

            for i = length(this.constraintLayers) : -1 : 1
                for j = 1:length(this.bodyLayers{i})
                    this.bodies{this.bodyLayers{i}(j)}.applyJacobiShock();
                end
                for iter = 1:this.iters
                    for j = 1 : length(this.constraintLayers{i})
                        this.constraintList{this.constraintLayers{i}(j)}.solve(this.hs);
                    end
                end
            end
            
            %{
            %this.draw();
            for i = length(this.bodyLayers):-1:1
                for j = 1:length(this.bodyLayers{i})
                    this.bodies{this.bodyLayers{i}(j)}.applyJacobiShock();
                end
                %this.draw();
                for iter = 1:50
                    for j = 1:length(this.bodyLayers{i})
                        body = this.bodies{this.bodyLayers{i}(j)};
                        for k = 1:length(body.shockParentConIndex)
                            collisionIndex = body.shockParentConIndex{k};
                            this.constraintList{collisionIndex}.solve(this.hs);
                        end
                    end
                end
            end
            %}
        %this.draw();
		end
		%%
		function solveConGSSP(this)
			%this.draw();
			%fprintf('substep %d\n',this.ks);
			for i = 1 : length(this.constraints)
				this.constraints{i}.clear();
            end
		    %this.collider.run();
			for iter = 0 : this.iters-1 % index 0
				%fprintf('  iter %d\n',iter);
				% Clear the Jacobi updates
				for i = 1 : length(this.bodies)
					this.bodies{i}.clearJacobi();
				end
				% Gauss-Seidel solve for non-collision constraints
				for j = 1 : length(this.constraints)
					this.constraints{j}.solve();
                end
            end

            collisions = this.collider.collisions;
            for i = 1 : length(collisions)
                collisions{i}.shockProp = true;
            end
		    for i = 1 : length(this.bodies)
                %this.draw();
                for iter = 1: 10
                    for j = 1 : length(this.bodies{i}.shockParentConIndex)
                        collisionIndex = this.bodies{i}.shockParentConIndex{j};
                        collisions{collisionIndex}.solve(this.hs);
                    end                  
                end
                %this.draw();
            end
                
            %this.draw();
            
            for i = 1 : length(collisions)
                collisions{i}.shockProp = true;
            end
		    for i = length(this.bodies) : -1 : 1
                this.bodies{i}.applyJacobiShock();
                %this.draw();
                for iter = 1: 10
                    for j = 1 : length(this.bodies{i}.shockParentConIndex)
                        collisionIndex = this.bodies{i}.shockParentConIndex{j};
                        collisions{collisionIndex}.solve(this.hs);
                    end  
                    %this.draw();
                end
                %this.draw();
            end

            %{
			for i = 1 : length(this.bodies)
                this.bodies{i}.x1_0 = this.bodies{i}.x;
			end
            %}
		end

		%%
		function solveConGS(this)
            %{
            %this.draw();
			%fprintf('substep %d\n',this.ks);
			for i = 1 : length(this.constraints)
				this.constraints{i}.clear();
            end
		    %this.collider.run();
            collisions = this.collider.collisions;

            for i = 1 : length(collisions)
                collisions{i}.shockProp = false;
            end
            %this.draw();
            for iter = 1 : this.iters
   				for j = 1 : length(this.constraints)
					this.constraints{j}.solve();
                end  

		        for i = 1 : length(this.bodies)
                    for j = 1 : length(this.bodies{i}.shockParentConIndex)
                        collisionIndex = this.bodies{i}.shockParentConIndex{j};
                        collisions{collisionIndex}.solve(this.hs)
                    end  
                end
            end
            %}

			for i = 1 : length(this.constraintList)
                this.constraintList{i}.shockProp = false;
            end
            
            for iter = 1 : this.iters
	            for i = 1 : length(this.constraintList)
                    this.constraintList{i}.solve(this.hs)  
                end
            end

        end

		%%
		function solveConGSO(this)
            %this.draw();
			%fprintf('substep %d\n',this.ks);
			for i = 1 : length(this.constraints)
				this.constraints{i}.clear();
            end
		    %this.collider.run();
            collisions = this.collider.collisions;
			for i = 1 : length(this.bodies)
				this.bodies{i}.clearJacobi();
			end
			for iter = 0 : this.iters-1 % index 0

				% Gauss-Seidel solve for non-collision constraints
				for j = 1 : length(this.constraints)
					this.constraints{j}.solve();
                end
                
                max_index = 1;
                max_d = 0;
				for j = 1 : length(collisions)
				    collisions{j}.update();
                    if abs(collisions{j}.d) > max_d
                        max_d = abs(collisions{j}.d);
                        max_index = j;
                    end
                end
                collisions{max_index}.solveNorPos();
                for i = 1 : length(this.bodies)
                        this.bodies{i}.applyJacobi();
                end
                %this.draw();

                %this.draw();
                %{
			    for j = length(collisions) : -1 : 1
                    collisions{j}.update();
				    collisions{j}.solveTanVel(this.k,this.ks,this.hs);
    			    for i = 1 : length(this.bodies)
				        this.bodies{i}.applyJacobi();
				        %fprintf('%e ',this.bodies{i}.x(end));
                    end
                    %this.draw()
				    %fprintf('%e\n',collisions{j}.C(1));
                end
                %this.draw()
                %}
            end
		end
		%%
		function computeEnergies(this)
			T = 0;
			V = 0;
			for i = 1 : length(this.bodies)
				[Ti,Vi] = this.bodies{i}.computeEnergies(this.k,this.ks,this.hs,this.grav);
				T = T + Ti;
				V = V + Vi;
			end
			for j = 1 : length(this.constraints)
				Vj = this.constraints{j}.computeEnergy();
				V = V + Vj;
			end
			this.tt(end+1) = this.t;
			this.TT(end+1) = T;
			this.VV(end+1) = V;
		end

		%%
		function draw(this)
			if this.drawHz == 0
				return;
			end
			if this.t == 0
				clf;
				xlabel('X');
				ylabel('Y');
				zlabel('Z');
				axis equal;
				if ~isempty(this.axis)
					axis(this.axis); %#ok<CPROP>
				end
				ax = gca;
				ax.Clipping = 'off';
				grid on;
				view(this.view); %#ok<CPROP>
			end
			if (floor(this.t*this.drawHz) > floor((this.t-this.h)*this.drawHz))
				cla;
				hold on;

				% Draw bodies
				nbodies = length(this.bodies);
				faces = cell(1,nbodies);
				verts = cell(1,nbodies);
				for i = 1 : nbodies
					[faces{i},verts{i}] = this.bodies{i}.draw();
				end

				% Draw constraints
				for i = 1 : length(this.constraints)
					this.constraints{i}.draw();
				end

				% Lighting
				l = light('Style','local','Position',[0 -50 100]);

				if this.ground.E(4,4) ~= 0 && this.ground.size > 0
					% Draw ground
					se3.drawAxis(this.ground.E);
					s = this.ground.size/2;
					V = this.ground.E(1:3,:)*[-s -s 0 1; s -s 0 1; s s 0 1; -s s 0 1]';
					F = [1 2 3 4];
					patch('Faces',F,'Vertices',V','FaceColor',[0.9 0.9 0.9]);

					% Shadow:
					%   p = l - ((d + dot(n,l))/dot(n,v - l))*(v - l),
					% where l is the light position, {n,d} is the plane, and v is the
					% vertex position.
					n = this.ground.E(1:3,3);
					d = -n'*this.ground.E(1:3,4);
					lpos = l.Position';
					for i = 1 : length(faces)
						F = faces{i};
						V = verts{i}';
						Vl = V - lpos;
						Vshadow = lpos - ((d + n'*lpos)./(n'*Vl)).*Vl + 1e-3*n;
						patch('Faces',F,'Vertices',Vshadow','EdgeColor','none','FaceColor',[0.2 0.2 0.2]);
					end
				end

				alpha(0.5);

				title(sprintf('t = %.4f',this.t));
				drawnow;

				if ~isempty(this.video)
					this.video.writeVideo(getframe(gcf));
				end
			end
		end

		%%
		function plotEnergy(this)
			this.VV = this.VV - this.VV(1);
			if this.plotH
				clf;
				hold on;
				plot(this.tt,this.TT,'-');
				plot(this.tt,this.VV,'-');
				plot(this.tt,this.TT+this.VV,'-');
				xlim([0,this.tEnd]);
				grid on;
				legend('T','V','H');
				xlabel('Time');
				ylabel('Energy');
			end
		end
	end

	%%
	methods (Static)
		%%
		function clearGlobal()
			global countB CM; %#ok<GVMIS>
			countB = 0;
			if ~usejava('jvm')
				CM = zeros(1,3);
			else
				CM = colormap('lines');
			end
		end

		%%
		function out = countB(incr)
			global countB; %#ok<GVMIS>
			if nargin < 1
				incr = 1;
			end
			countB = countB + incr;
			out = countB;
		end

		%%
		function test(modelID)
			if nargin < 1
				modelID = 0;
			end
			model = createTestModels(modelID);
			model.init();
			model.simulate();
			model.plotEnergy();
		end
	end
end
