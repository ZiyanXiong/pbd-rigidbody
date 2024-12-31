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
                %this.solveConGlobal();
                this.solveConGPQP();
                %this.solveCon();
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
                        this.collider.collisions{j}.solveCollisionNor(-Inf);
                        this.collider.collisions{j}.solveCollisionTan();
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
                        this.collider.collisions{j}.solveCollisionNor(0);
                        this.collider.collisions{j}.solveCollisionTan();
                    end
                end
            end

            for i = 1 : length(this.bodies)
                this.bodies{i}.integrateStates();
            end
            
        end

		%%
        function solveConGPQP(this)
            this.collider.run();

            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    this.collider.collisions{j}.initConstraints(this.h, this.hs);
                end
            end
            this.draw();

			this.stepBDF1();

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

            output = this.GPQP(n);
            lambdas = output.lambdas;
            %lambdas = solver.Gauss_Sidiel(A, b, mu);
            
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

			this.t = this.t + this.h;
            for i = 1 : length(this.bodies)
                this.bodies{i}.integrateStates();
            end

        end

		%%
        function solveConGlobal(this)
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

                
                mu = this.bodies{1}.mu;
                itermax = 1000;                
                solver = ConstraintSolver(itermax,1e-6);
                lambdas = solver.Cone_GPQP(A, b, mu);
                %lambdas = solver.Gauss_Sidiel(A, b, mu);
                
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

        function output = GPQP(this, n)
            tol = 1e-8;
            eps = 1e-8;
            iterMax = this.substeps;   
            CGiterMax = 200;
            rs = zeros(iterMax,1);
            collisions = [];
            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    collisions(end+1) = j;
                end
            end


            fPrev = 0;
            gPrev = zeros(n,1);
            g = zeros(n,1);
            lambda = zeros(n,1);
            
            for i = 1:length(this.bodies)
                this.bodies{i}.LTx = zeros(6,1);
            end
            for i = collisions
                this.collider.collisions{i}.compute_LTlambda();
            end
            for i = collisions
                this.collider.collisions{i}.compute_LLTx();
                this.collider.collisions{i}.g = this.collider.collisions{i}.Ax - this.collider.collisions{i}.b;
                fPrev = fPrev + this.collider.collisions{i}.lambda' * ( 0.5 * this.collider.collisions{i}.Ax -  this.collider.collisions{i}.b);
                gPrev(this.collider.collisions{i}.mIndces) = this.collider.collisions{i}.g;
            end

            CGiterVec = [];
            
            for iter = 1:iterMax
                if(iter == 3)
                    disp(iter);
                end

                % Compute Cauchy Point
                tList = Inf(n, 1);
                for i = 1:length(this.bodies)
                    this.bodies{i}.LTx = zeros(6,1);
                end
                for i = collisions
                    this.collider.collisions{i}.compute_LTlambda();
                end
                for i = collisions
                    this.collider.collisions{i}.compute_LLTx();
                    this.collider.collisions{i}.g = this.collider.collisions{i}.Ax - this.collider.collisions{i}.b;
                    tList(this.collider.collisions{i}.mIndces) = this.collider.collisions{i}.compute_tbar();
                end

                tUniqueList = unique(tList,'sorted');
                if(tUniqueList(end) ~= Inf)
                    tUniqueList(end + 1) = Inf;
                end
                tc = 0;
                for tIndex = 1: size(tUniqueList,1)
                    for i = collisions
                        this.collider.collisions{i}.compute_p(tc);
                        this.collider.collisions{i}.compute_lambdac(tc);
                    end

                    fPrime = 0;
                    fPrimePrime = 0;

                    for i = 1:length(this.bodies)
                        this.bodies{i}.LTx = zeros(6,1);
                    end
                    for i = collisions
                        this.collider.collisions{i}.compute_LTp();
                    end
                    for i = collisions
                        this.collider.collisions{i}.compute_LLTx();
                        fPrime = fPrime - this.collider.collisions{i}.b' * this.collider.collisions{i}.p + ...
                                    this.collider.collisions{i}.lambdac' * this.collider.collisions{i}.Ax;
                        fPrimePrime = fPrimePrime + this.collider.collisions{i}.p' * this.collider.collisions{i}.Ax;
                    end

                    deltaTStar = - fPrime / fPrimePrime;
                    if(fPrime >0)
                        break;
                    elseif(deltaTStar >=0 && deltaTStar < tUniqueList(tIndex) - tc)
                        tc = tc + deltaTStar;
                        break;
                    end
                    tc = tUniqueList(tIndex);
                end

                for i = collisions
                    this.collider.collisions{i}.compute_lambdac(tc);
                    this.collider.collisions{i}.compute_lambdad(tc);
                    this.collider.collisions{i}.lambda = this.collider.collisions{i}.lambdac;
                end

                % PCG
                % Compute b_cg
                for i = collisions
                     this.collider.collisions{i}.compute_degenerate_J1I_J2I_b();
                end

                % Init CG 
                for i = 1:length(this.bodies)
                    this.bodies{i}.LTx = zeros(6,1);
                end
                for i = collisions
                    this.collider.collisions{i}.compute_degenerate_LTlambda();
                end
                for i = collisions
                    this.collider.collisions{i}.compute_degenerate_LLTx();
                    this.collider.collisions{i}.r_cg = this.collider.collisions{i}.Ax - this.collider.collisions{i}.b_cg;
                    this.collider.collisions{i}.g_cg = this.collider.collisions{i}.Minv_cg .* this.collider.collisions{i}.r_cg;
                    this.collider.collisions{i}.d_cg = - this.collider.collisions{i}.g_cg;
                end
                % CG iterations 
                for CGiter = 1: CGiterMax
                    r = 0;
                    for i = collisions
                        M = 1 ./ this.collider.collisions{i}.Minv_cg;
                        r = r + this.collider.collisions{i}.r_cg(this.collider.collisions{i}.freeIndex)' * diag(M(this.collider.collisions{i}.freeIndex)) ...
                            * this.collider.collisions{i}.r_cg(this.collider.collisions{i}.freeIndex);
                    end
                    if(r < tol)
                        break;
                    end

                    numerator = 0;
                    denominator = 0;

                    for i = 1:length(this.bodies)
                        this.bodies{i}.LTx = zeros(6,1);
                    end
                    for i = collisions
                        this.collider.collisions{i}.compute_LTd_cg();
                    end
                    for i = collisions
                        this.collider.collisions{i}.compute_degenerate_LLTx();
                        numerator = numerator + this.collider.collisions{i}.r_cg(this.collider.collisions{i}.freeIndex)' ...
                                    * this.collider.collisions{i}.g_cg(this.collider.collisions{i}.freeIndex);
                        denominator = denominator + this.collider.collisions{i}.d_cg(this.collider.collisions{i}.freeIndex)' ...
                                        * this.collider.collisions{i}.Ax(this.collider.collisions{i}.freeIndex);
                    end
                    alpha = numerator / denominator;
                
                    feasible = true;
                    for i = collisions
                        feasible = feasible & this.collider.collisions{i}.update_cg(alpha);
                    end
                    if(~feasible)
                        break;
                    end

                    numerator = 0;
                    denominator = 0;
                    for i = collisions
                        denominator = denominator + this.collider.collisions{i}.r_cg(this.collider.collisions{i}.freeIndex)' ...
                                        * this.collider.collisions{i}.g_cg(this.collider.collisions{i}.freeIndex);
                        this.collider.collisions{i}.r_cg = this.collider.collisions{i}.r_cg + alpha * this.collider.collisions{i}.Ax;
                        this.collider.collisions{i}.g_cg = this.collider.collisions{i}.Minv_cg .* this.collider.collisions{i}.r_cg;
                        numerator = numerator + this.collider.collisions{i}.r_cg(this.collider.collisions{i}.freeIndex)' ...
                                    * this.collider.collisions{i}.g_cg(this.collider.collisions{i}.freeIndex);
                    end
                    beta = numerator / denominator;
                    for i = collisions
                        this.collider.collisions{i}.d_cg = -this.collider.collisions{i}.g_cg + beta * this.collider.collisions{i}.d_cg;
                    end
                end

                CGiterVec = [CGiterVec, CGiter];
            
                % Project to feasible region
                for i = collisions
                    this.collider.collisions{i}.project();
                end

                % Test if results satisfy the KKT conditions
                f = 0;
                for i = 1:length(this.bodies)
                    this.bodies{i}.LTx = zeros(6,1);
                end
                for i = collisions
                    this.collider.collisions{i}.compute_LTlambda();
                end
                for i = collisions
                    this.collider.collisions{i}.compute_LLTx();
                    this.collider.collisions{i}.g = this.collider.collisions{i}.Ax - this.collider.collisions{i}.b;
                    f = f + this.collider.collisions{i}.lambda' * ( 0.5 * this.collider.collisions{i}.Ax -  this.collider.collisions{i}.b);
                    g(this.collider.collisions{i}.mIndces) = this.collider.collisions{i}.g;
                    lambda(this.collider.collisions{i}.mIndces) = this.collider.collisions{i}.lambda;
                end

	            if norm(g) < eps
		            break;
	            end
                if norm(g - gPrev) < eps
                    break;
                end
                if norm(f-fPrev) < eps
                    break;
                end
                fPrev = f;
                gPrev = g;
                rs(iter) = norm(g);
            end
            
            
            output.iterations = iter;
            output.lambdas = lambda;
            output.cgiterations = CGiterVec;
            output.rs = rs;
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
