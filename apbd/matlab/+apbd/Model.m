classdef Model < handle
	%Model Class to hold model data

	properties
		bodies % list of bodies
		constraints %list of constraints
		collider % collision handler
		grav % gravity
		ground % ground transform, with Z up
        biasCoefficient
		
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
        solverType % 1: PhysX 2:PhysX+OneWayShockPropagation 3:PhysX+TwoWayShockPropagation
        savedBodyStatesPath
        
	end

	methods
		function this = Model()
			apbd.Model.clearGlobal();
			
			% Default values
			this.bodies = {};
			this.constraints = {};
			this.grav = [0 0 -980]';
			this.ground.E = zeros(4);

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
            this.solverType = 1;
		end

		%%
		function init(this)
            this.collider = apbd.Collider(this);
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
            if ~isempty(this.savedBodyStatesPath)
                this.saveBodyStates();
            else
                this.draw();
            end
		end

		%%
		function simulate(this)
			while this.k < this.steps
				this.ks = 0;
                if this.k == 4
                    %fprintf("Pause.");
                end
                %this.solveConTGS();
                this.solveCon();
				this.k = this.k + 1;
				this.computeEnergies();
				%this.draw();
                if ~isempty(this.savedBodyStatesPath)
                    this.saveBodyStates();
                else
                    this.draw();
                end
			end
        end

		%%
		function stepBDF1(this)
			for i = 1 : length(this.bodies)
				this.bodies{i}.stepBDF1(this.h,this.grav);
			end
        end

		%%
		function solveConTGS(this)
			this.stepBDF1();
            this.collider.run();
            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    this.collider.collisions{j}.initConstraints(this.h, this.hs);
                end
            end
            this.draw();
            
            if(this.solverType > 1)
                % Shock Propagation
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        for k = 1:5
                            this.collider.collisions{j}.solveCollisionNor(-Inf, true);
                            %this.collider.collisions{j}.solveCollisionTan(true);
                        end
                    end
                end
                if(this.solverType > 2)
                    for i = length(this.collider.activeCollisions) : -1 : 1
                        for j = this.collider.activeCollisions{i}
                            for k = 1:25
                                this.collider.collisions{j}.solveCollisionNor(-Inf, true);
                                %this.collider.collisions{j}.solveCollisionTan(true);
                            end
                        end
                        for j = this.collider.activeCollisions{i}
                            this.collider.collisions{j}.applyLambdaSP();
                        end
                    end
                end

                for i = 1 : length(this.bodies)
                    %this.bodies{i}.updateStates(this.hs);
                end 
            end

            while this.ks < this.substeps
			    %this.draw();
			    %fprintf('substep %d\n',this.ks);
                for i = 1 : length(this.constraints)
				    this.constraints{i}.clear();
                end

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

                
                % Gauss-Seidal 
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        this.collider.collisions{j}.solveCollisionNor(-Inf, false);
                        this.collider.collisions{j}.solveCollisionTan(false);
                    end
                end
                
                for i = 1 : length(this.bodies)
                    this.bodies{i}.updateStates(this.hs);
                end
                %this.draw();
				this.t = this.t + this.hs;
				this.ks = this.ks + 1;
            end

            for iter = 1: 1
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        this.collider.collisions{j}.solveCollisionNor(0, false);
                        this.collider.collisions{j}.solveCollisionTan(false);
                    end
                end
            end

            for i = 1 : length(this.bodies)
                this.bodies{i}.integrateStates();
            end
            
        end


		%%
        function solveCon(this)
            this.collider.run();

            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    this.collider.collisions{j}.initConstraints(this.h, this.hs);
                end
            end
            this.draw();

			this.stepBDF1();
            for iter = 1 : this.iters
			    %this.draw();
			    %fprintf('substep %d\n',this.ks);
                for i = 1 : length(this.constraints)
				    this.constraints{i}.clear();
                end

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

                n = 1;
                ci = 1;
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        this.collider.collisions{j}.index = ci;
                        this.collider.collisions{j}.mIndces = n : n - 1 + this.collider.collisions{j}.contactNum * 3;
                        n = n + this.collider.collisions{j}.contactNum * 3;
                        this.collider.collisions{j}.computeJ_b();
                        ci = ci + 1;
                    end
                end

                n = n-1;
                A = zeros(n,n);
                Asp = zeros(n,n);
                b = zeros(n,1);
                blocks = {};
                % GP
                for i = 1 : length(this.bodies)
                    for j = 1 : length(this.bodies{i}.collisions)
                        cj = this.bodies{i}.collisions(j);
                        this.collider.collisions{cj}.nextColl = [];
                        for k = j + 1 : length(this.bodies{i}.collisions)
                            ck = this.bodies{i}.collisions(k);
                            this.collider.collisions{cj}.nextColl(end+1) = ck;
                            rows = this.collider.collisions{cj}.mIndces;
                            cols = this.collider.collisions{ck}.mIndces;
                            if(this.collider.collisions{cj}.layer == this.collider.collisions{ck}.layer)
                                A(rows,cols) = this.collider.collisions{cj}.J1I * this.collider.collisions{ck}.J1I';
                                %Asp(cols,rows) = A(rows,cols)';
                            else
                                A(rows,cols) = this.collider.collisions{cj}.J1I * this.collider.collisions{ck}.J2I';
                            end
                            Asp(rows,cols) = A(rows,cols);
                            A(cols,rows) = A(rows,cols)';
                        end
                    end
                end

                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        inds = this.collider.collisions{j}.mIndces;
                        A(inds,inds) = this.collider.collisions{j}.J1I * this.collider.collisions{j}.J1I' + this.collider.collisions{j}.J2I * this.collider.collisions{j}.J2I';
                        Asp(inds,inds) = this.collider.collisions{j}.J1I * this.collider.collisions{j}.J1I';
                        %Asp(inds,inds) = A(inds,inds);
                        b(inds) = this.collider.collisions{j}.b;
                        blocks{end+1} = inds;
                    end
                end

                L = zeros(n,6*length(this.bodies));
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        rows = this.collider.collisions{j}.mIndces;
                        if(this.collider.collisions{j}.ground)
                            cols =  (this.collider.collisions{j}.body1.index - 1) * 6 + 1 : this.collider.collisions{j}.body1.index * 6;
                            L(rows,cols) = this.collider.collisions{j}.J1I;
                        else
                            cols =  (this.collider.collisions{j}.body1.index - 1) * 6 + 1 : this.collider.collisions{j}.body1.index * 6;
                            L(rows,cols) = this.collider.collisions{j}.J1I;
                            cols =  (this.collider.collisions{j}.body2.index - 1) * 6 + 1 : this.collider.collisions{j}.body2.index * 6;
                            L(rows,cols) = this.collider.collisions{j}.J2I;
                        end
                    end
                end
                A = L*L';
                Lsp = L;
                blocks = {};
                for i = 1 : length(this.collider.activeCollisions)
                    inds = [];
                    for j = this.collider.activeCollisions{i}
                        rows = this.collider.collisions{j}.mIndces;
                        inds = [inds, rows];
                        if(~this.collider.collisions{j}.ground && mod(this.collider.collisions{j}.layer,2) == 1)
                            cols =  (this.collider.collisions{j}.body2.index - 1) * 6 + 1 : this.collider.collisions{j}.body2.index * 6;
                            Lsp(rows,cols) = 0;
                        end
                    end
                    if(~isempty(inds)&& mod(this.collider.collisions{j}.layer,2) == 0)
                        blocks{end} = [blocks{end}, inds];
                    else
                        blocks{end+1} = inds;
                    end
                end
                Asp = L * Lsp';
                %{
                Aspinv = Asp;
                for i = 2:5
                    Aspinv(blocks{i},blocks{i})= A(blocks{i},blocks{i}) - A(blocks{i-1},blocks{i})' * pinv(Aspinv(blocks{i-1},blocks{i-1})) *  A(blocks{i-1},blocks{i});
                end
                %}


                mu = this.bodies{1}.mu;
                itermax = 1000;
                
                clf;
                solver = ConstraintSolver(itermax,1e-6);
                lambdas = solver.Gauss_Sidiel(A, b, mu);
                solver.draw('Gauss-Seidel');
                
                load('SP_allproj.mat');
                solver = ConstraintSolver(itermax,1e-6);
                lambdas = solver.Gauss_Sidiel(A, b, mu, lambdas);
                solver.draw('Gauss-Seidel with SP warm start (projection on noraml and tangent dirctions)');

                load('SP_norproj.mat');
                solver = ConstraintSolver(itermax,1e-6);
                lambdafull = zeros(length(b),1);
                lambdafull(1:3:end) = lambdas;
                lambdas = solver.Gauss_Sidiel(A, b, mu, lambdafull);
                solver.draw('Gauss-Seidel with SP warm start (projection on noraml direction)');

                load('SP_noproj.mat');
                solver = ConstraintSolver(itermax,1e-6);
                lambdas = solver.Gauss_Sidiel(A, b, mu, lambdas);
                solver.draw('Gauss-Seidel with SP warm start (No projection)');
                %}

                %lambdas = solver.Chel(A, b, mu);
                %solver.draw('Chel');
                %lambdas = solver.Shock_Propagation(A, b, Asp, blocks, mu);
                %solver.draw('Shock-Propagation');
                %save('SP_norproj.mat', "lambdas");
                %lambdas = solver.Preconditioned_Conjugate_Gradient(A, b, Asp, mu);
                %solver.draw('Preconditioned-Conjugate-Gradient');
                %lambdas = solver.Conjugate_Gradient(A, b, mu);
                %solver.draw('Conjugate-Gradient');
                %lambdas = solver.Mix_GS_CG(A, b, mu);
                %solver.draw('Mixing Gauss-Seidel and Conjugate-Gradient');
                legend
                xlabel('Iteration number') 
                ylabel('Error') 
                title('Gauss-Seidel solver for arch scene with different warm start')
                %save('arch.mat',"A","b");
                
                solver = ConstraintSolver(itermax,1e-6);
                lambdas = solver.Gauss_Sidiel(A, b, mu);
                %lambdas = solver.Shock_Propagation(A, b, Asp, blocks, mu);
                lambdas = pinv(A)*b;
                %lambdas(1:end) = 0;
                
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        l = this.collider.collisions{j}.contactNum*3 - 1;
                        start = this.collider.collisions{j}.mIndces(1);
                        lambdai = lambdas(start:start+l);
                        for k = 1: this.collider.collisions{j}.contactNum
                            this.collider.collisions{j}.constraints{k}.applyLambda(lambdai(3*(k-1) + 1: 3*k));
                        end
                    end
                end

                for i = 1 : length(this.bodies)
                    this.bodies{i}.updateStatesDirect(this.h);
                end
            end

			this.t = this.t + this.h;
            for i = 1 : length(this.bodies)
                this.bodies{i}.integrateStates();
            end

        end
        %%
        function dfdp = backward(this)
            x_target = [0 0 0 1 0 0 0]';
            stateLength = 6;
            bodyNum = length(this.bodies);
            conNum = length(this.collider.activeCollisions);
            n = this.steps * bodyNum * stateLength;
            y = zeros(n,1);
            H = eye(n);
            xs = [];
            for i  = 1:length(this.bodies{1}.xs_bar)
                for j = 1:bodyNum
                    xs = [xs;this.bodies{j}.xs_bar{i}];
                end
            end
            y(this.steps*stateLength:end) = xs(this.steps*stateLength:end) - x_target;  

            for i = 1 : this.steps
                if i == 2
                    for p = 1 : bodyNum
                        this.bodies{p}.x = xs(ind.getInd(i-1,this.iters,conNum));
                        this.bodies{p}.x0 = this.bodies{p}.xInit;
                        index = this.bodies{p}.index;
		                [dUncondu,~] = this.bodies{p}.dUnconudPrevu(this.hs);
                        H(ind.getInd(i,0,0,index),ind.getInd(i-1,this.iters,conNum,index)) = -dUncondu;
                    end
                end

                if i > 2
                    for p = 1 : bodyNum
                        this.bodies{p}.x = xs(ind.getInd(i-1,this.iters,conNum));
                        this.bodies{p}.x0 =  xs(ind.getInd(i-2,this.iters,conNum));
                        index = this.bodies{p}.index;
		                [dUncondu,dUncondu0] = this.bodies{p}.dUnconudPrevu(this.hs);
                        H(ind.getInd(i,0,0,index),ind.getInd(i-1,this.iters,conNum,index)) = -dUncondu;
                        H(ind.getInd(i,0,0,index),ind.getInd(i-2,this.iters,conNum,index)) = -dUncondu0;
                    end
                end

                for j = 1 : this.iters
                    for p = 1 : conNum
                        if length(this.constraintList{p}.bodies) == 1
                            index = this.constraintList{p}.body.index;
                            this.constraintList{p}.body.x = xs(ind.getInd(i,j,p-1,index));
                            this.constraintList{p}.body.x1_0 = xs(ind.getInd(i,0,0,index));
		                    [dnextdx,dnextdx0] = this.constraintList{p}.dnextdx();
                            H(ind.getInd(i,j,p,index),ind.getInd(i,j,p-1,index)) = H(ind.getInd(i,j,p,index),ind.getInd(i,j,p-1,index)) - dnextdx;
                            H(ind.getInd(i,j,p,index),ind.getInd(i,0,0,index)) = H(ind.getInd(i,j,p,index),ind.getInd(i,0,-1,index)) - dnextdx0;
                            for uindex = 1:bodyNum
                                if uindex ~= index
                                    H(ind.getInd(i,j,p,uindex),ind.getInd(i,j,p-1,uindex)) = -eye(7);
                                end
                            end
                        end
                        if length(this.constraintList{p}.bodies) == 2
                            index1 = this.constraintList{p}.body1.index;
                            index2 = this.constraintList{p}.body2.index;
                            this.constraintList{p}.body1.x = xs(ind.getInd(i,j,p-1,index1));
                            this.constraintList{p}.body1.x1_0 = xs(ind.getInd(i,0,0,index1));
                            this.constraintList{p}.body2.x = xs(ind.getInd(i,j,p-1,index2));
                            this.constraintList{p}.body2.x1_0 = xs(ind.getInd(i,0,0,index2));
                            this.constraintList{p}.setNormal();
                            this.constraintList{p}.testGradient();
		                    [dnext1dx1, dnext1dx2, dnext1dx10, dnext1dx20, dnext2dx1, dnext2dx2, dnext2dx10, dnext2dx20] = this.constraintList{p}.dnextdx();
                            H(ind.getInd(i,j,p,index1),ind.getInd(i,j,p-1,index1)) = H(ind.getInd(i,j,p,index1),ind.getInd(i,j,p-1,index1)) - dnext1dx1;
                            H(ind.getInd(i,j,p,index1),ind.getInd(i,0,0,index1)) = H(ind.getInd(i,j,p,index1),ind.getInd(i,0,-1,index1)) - dnext1dx10;
                            H(ind.getInd(i,j,p,index1),ind.getInd(i,j,p-1,index2)) = H(ind.getInd(i,j,p,index1),ind.getInd(i,j,p-1,index2)) - dnext1dx2;
                            H(ind.getInd(i,j,p,index1),ind.getInd(i,0,0,index2)) = H(ind.getInd(i,j,p,index1),ind.getInd(i,0,-1,index2)) - dnext1dx20;

                            H(ind.getInd(i,j,p,index2),ind.getInd(i,j,p-1,index1)) = H(ind.getInd(i,j,p,index2),ind.getInd(i,j,p-1,index1)) - dnext2dx1;
                            H(ind.getInd(i,j,p,index2),ind.getInd(i,0,0,index1)) = H(ind.getInd(i,j,p,index2),ind.getInd(i,0,-1,index1)) - dnext2dx10;
                            H(ind.getInd(i,j,p,index2),ind.getInd(i,j,p-1,index2)) = H(ind.getInd(i,j,p,index2),ind.getInd(i,j,p-1,index2)) - dnext2dx2;
                            H(ind.getInd(i,j,p,index2),ind.getInd(i,0,0,index2)) = H(ind.getInd(i,j,p,index2),ind.getInd(i,0,-1,index2)) - dnext2dx20;
                        end
                    end
                end
            end
            z = H'\ y;
            %dfdp = z(ind.getInd(1,0,0)) - z(ind.getInd(2,0,0));
            this.bodies{1}.x = this.bodies{1}.xInit;
            this.bodies{1}.x0 = this.bodies{1}.xInit;
            [dUncondu,dUncondu0] = this.bodies{1}.dUnconudPrevu(this.hs);
            dfdp = (dUncondu + dUncondu0)' * z(ind.getInd(1,0,0,1));
            %{
            this.bodies{1}.x = xs(ind.getInd(1,this.iters,conNum));
            this.bodies{1}.x0 = this.bodies{1}.xInit;
            [~,dUncondu0] = this.bodies{1}.dUnconudPrevu(this.hs);
            dfdp = dfdp + dUncondu0' * z(ind.getInd(2,0,0,1));
            %}
        end

        %%
        function H = getH_FD(this)

        end
        %%
        function [H_prev, H] = solveTGSGradient(this)
            stateLength = 6;
            constraintNum = 4;
            n = (2 + this.substeps * constraintNum);
            H_prev = zeros(n * stateLength,n * stateLength);
            H = ones(n * stateLength,n * stateLength);
            H_prev(1:stateLength, (n-1) * stateLength:end) = - ones(stateLength);
            for iter = 1 : this.substeps
                % Gauss-Seidal 
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        H_sub = this.collider.collisions{j}.getGradient(-Inf, false);
                        H = H + H_sub;
                    end
                end
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

				% Draw collisions
		        for i = 1 : length(this.collider.activeCollisions)
                    for collisions = this.collider.activeCollisions
                        for j = collisions{1}
                            this.collider.collisions{j}.draw();
                        end
                    end
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

        function saveBodyStates(this)
            % Check if the file exists
            if exist(this.savedBodyStatesPath, 'file') == 0
                % If the file does not exist, create it
                fid = fopen(this.savedBodyStatesPath, 'w');
                fclose(fid);
            end
            
            if (floor(this.t*this.drawHz) > floor((this.t-this.h)*this.drawHz))
                % Append new lines to the file
                fid = fopen(this.savedBodyStatesPath, 'a'); % 'a' mode for appending
                % Iterate through each quaternion and position
                for i = 1 : length(this.bodies)
                    fprintf(fid, '%f %f %f %f %f %f %f ', this.bodies{i}.x);
                end
                fprintf(fid, '\n');
                fclose(fid);
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
