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
        solverType % 1: GS 2: 2PSP 3: TGS
        savedBodyStatesPath
        iterVec
        rVec
        modelID
        resultFolder
        useContactCaching
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
            this.useContactCaching = false;
			this.Hexpected = zeros(1,2);

			% Drawing etc.
			this.name = 'Default';
			this.ground.size = 10;
			this.drawHz = 15;
			this.view = 3;
			this.axis = [];
			this.video = [];
            this.solverType = 1;
            this.iterVec = [];
            this.rVec = [];
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
                if this.k == 2
                    fprintf("Pause.");
                end
                if(this.solverType == 3)
                    this.solveConTGS();
                else
                    this.solveConGlobal();
                end
                %this.solveConGPQP();
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
            f = zeros(3,1);
			for i = 1 : length(this.bodies)
                if(this.modelID == 11 && i==8 && this.k<20)
                    this.bodies{i}.stepBDF1(this.h,this.grav,[0 -800000 0]');
                else
				    this.bodies{i}.stepBDF1(this.h,this.grav,f);
                end
			end
        end

		%%
		function solveConTGS(this)
            this.collider.run();
            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    this.collider.collisions{j}.initConstraints(this.h, this.hs);
                end
            end
            this.draw();

            if(this.modelID == 9)
                groundVelocity = zeros(6,1);
                groundVelocity(4) = 8 * sin(this.k / 15 *pi);
                groundVelocity(6) = 50 * sin(this.k / 12 *pi);
                this.bodies{1}.setInitVelocity(groundVelocity);
            end

			this.stepBDF1();
            while this.ks < this.substeps
                %this.bodies{i}.stepBDF1(this.hs,this.grav,f);
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
                    
                    for j = this.collider.activeCollisions{i}
                        %this.collider.collisions{j}.solveCollisionTan();
                    end
                    
                end
                
                for i = 1 : length(this.collider.activeCollisions)
                    for j = this.collider.activeCollisions{i}
                        %this.collider.collisions{j}.solveCollisionTan();
                    end
                end

                for i = 1 : length(this.bodies)
                    this.bodies{i}.updateStates(this.hs);
                end
                %this.draw();
				this.t = this.t + this.hs;
				this.ks = this.ks + 1;
            end

            for iter = 1: 0
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

            n = 1;
            ci = 1;
            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    this.collider.collisions{j}.index = ci;
                    this.collider.collisions{j}.mIndces = n : n - 1 + this.collider.collisions{j}.contactNum * 3;
                    n = n + this.collider.collisions{j}.contactNum * 3;
                    this.collider.collisions{j}.compute_b();
                    ci = ci + 1;
                end
            end
            rs = zeros(n-1,1);
            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    rs(this.collider.collisions{j}.mIndces) = this.collider.collisions{j}.b;
                end
            end

            this.iterVec(end+1) = this.substeps;
            this.rVec(end+1) = norm(rs(rs>0));
            fid = fopen(fullfile(this.resultFolder, sprintf('Body_States_TGS_%d.txt',this.substeps)), 'a+');
            if(this.k==0)
                fprintf(fid, '#Body number: %d\n', length(this.bodies));
                fprintf(fid, '#Step number: %d\n', this.steps);
            end
            for i = 1:length(this.bodies)
                fprintf(fid, '%f %f %f %f %f %f %f ', this.bodies{i}.x);
            end
            fprintf(fid, '\n');
            fclose(fid);
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

            if(this.modelID == 9)
                groundVelocity = zeros(6,1);
                groundVelocity(4) = 5 * sin(this.k / 15 *pi);
                groundVelocity(6) = 65* sin(this.k / 12 *pi);
                this.bodies{1}.setInitVelocity(groundVelocity);
            end
            
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
                        this.collider.collisions{j}.compute_d();
                        ci = ci + 1;
                    end
                end

                n = n-1;
                A = zeros(n,n);
                Asp = zeros(n,n);
                b = zeros(n,1);
                d = zeros(n,1);
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
                        d(inds) = this.collider.collisions{j}.d;
                    end
                end

                [~, ind] = sort(cellfun(@(body) body.layer, this.bodies));
                for i = 1:length(this.bodies)
                    this.bodies{ind(i)}.colIndices = (i - 1)*6+1 : i*6;
                end

                L = zeros(n,6*length(this.bodies));
                L_sp = zeros(n,6*length(this.bodies));
                blocks = {};
                for i = 1 : length(this.collider.activeCollisions)
                    block = [];
                    for j = this.collider.activeCollisions{i}
                        rows = this.collider.collisions{j}.mIndces;
                        block = [block rows];
                        if(this.collider.collisions{j}.ground)
                            cols =  this.collider.collisions{j}.body1.colIndices;
                            L(rows,cols) = this.collider.collisions{j}.J1I;
                            L_sp(rows,cols) = this.collider.collisions{j}.J1I;
                        else
                            cols =  this.collider.collisions{j}.body1.colIndices;
                            L(rows,cols) = this.collider.collisions{j}.J1I;
                            L_sp(rows,cols) = this.collider.collisions{j}.J1I;
                            cols =  this.collider.collisions{j}.body2.colIndices;
                            L(rows,cols) = this.collider.collisions{j}.J2I;
                            if(this.collider.collisions{j}.body1.layer == this.collider.collisions{j}.body2.layer)
                                L_sp(rows,cols) = this.collider.collisions{j}.J2I;
                            end
                        end
                    end
                    if(this.collider.collisions{j}.body1.layer == this.collider.collisions{j}.body2.layer)
                        if(isempty(blocks))
                            blocks{end+1} = block;
                        else
                            blocks{end} = [blocks{end} block];
                        end
                    else
                        blocks{end+1} = block;
                    end
                end
                A = L*L';
                Asp = L * L_sp';
                
                mu = this.bodies{1}.mu;
                itermax = 150;                
                solver = ConstraintSolver(itermax,1e-6);
                %lambdas = solver.Cone_GPQP(A, b, mu);
                %lambdas = solver.SOCP(L, b, mu);
                %lambdas = solver.Gauss_Sidiel(A, b, mu);
                %lambdas = solver.Shock_Propagation_re(A, Asp, b, blocks, mu);
                %lambdas = solver.Shock_Propagation(A, b, Asp, mu);
                if(this.solverType ==  1)
                    lambdas = solver.Temporal_Gauss_Sidiel(A, b, mu, d, this.substeps);
                else
                    lambdas = solver.Shock_Propagation_lbl(A, Asp, b, d, blocks, mu);
                end
                this.iterVec(end+1) = solver.itercount;
                this.rVec(end+1) = solver.rs(this.iterVec(end));
                %{
                if(this.solverType == 2 && this.k == 1)
                    %clf;
                    solver.Gauss_Sidiel(A, b, mu);
                    %solver.draw('Gauss-Seidel, t = 1e-2');
                    rs_gs = solver.rs;
                    
                    %{
                    %save("gs_rs_1e-4.mat", "rs_gs_1e_4");
                    load("gs_rs_1e-3.mat");
                    semilogy(1:size(rs_gs_1e_3,1), rs_gs_1e_3, 'DisplayName','Gauss-Seidel, t = 1e-3','linewidth',2);
                    hold on;
        
                    load("gs_rs_1e-4.mat");
                    semilogy(1:size(rs_gs_1e_4,1), rs_gs_1e_4, 'DisplayName','Gauss-Seidel, t = 1e-4','linewidth',2);
                    hold on;
                    %}
        
                    %lambdas = solver.Staggered(A, b, mu);
                    %solver.draw('Staggered Projections');
                    %lambdas = solver.Shock_Propagation(A, b, Asp, mu);
                    solver.Shock_Propagation_lbl(A, Asp, b, blocks, mu);
                    %lambdas = solver.Shock_Propagation_Mix(L,L_sp, b, blocks, mu);
                    %lambdas = solver.Shock_Propagation_re(A,Asp, b, blocks, mu);
                    %solver.draw('2 Way Shock Propagation, t = 1e-2');
                    rs_2psp = solver.rs;
                    %lambdas = solver.Cone_GPQP(A, b, mu);
                    %solver.draw('Gradient Projection QP');
                    %legend
                    %xlabel('Iteration number') 
                    %ylabel('Residual') 
                    %title('Convergence plot for different solver')
                    save(strcat(this.resultFolder,sprintf("residual_per_iteration\\2PSP_%d.mat",this.k)), "rs_gs", "rs_2psp");
                end
                %}
                %{
                if(this.k == 30)
                    E_List = {};
                    for i = 1:length(this.bodies)
                        E_List{end+1} = this.bodies{i}.computeTransform();
                    end
                    save("Initial_Transform.mat","E_List");
                end
                %}
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
            if(this.solverType == 1)
                 fid = fopen(fullfile(this.resultFolder, sprintf('Body_States_TGS_%d.txt',this.substeps)), 'a+');
                if(this.k==0)
                    fprintf(fid, '#Body number: %d\n', length(this.bodies));
                    fprintf(fid, '#Step number: %d\n', this.steps);
                end
                for i = 1:length(this.bodies)
                    fprintf(fid, '%f %f %f %f %f %f %f ', this.bodies{i}.x);
                end
                fprintf(fid, '\n');
                fclose(fid);
            else
                fid = fopen(fullfile(this.resultFolder, 'Body_States_2PSP.txt'), 'a+');
                if(this.k==0)
                    fprintf(fid, '#Body number: %d\n', length(this.bodies));
                    fprintf(fid, '#Step number: %d\n', this.steps);
                end
                for i = 1:length(this.bodies)
                    fprintf(fid, '%f %f %f %f %f %f %f ', this.bodies{i}.x);
                end
                fprintf(fid, '\n');
                fclose(fid);
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

        function save_states(this,A,b)
            bodyList = {};
            for i = 1 : length(this.bodies)
                body.m = this.bodies{i}.Mp;
                body.I = this.bodies{i}.Mr;
                bodyList{end+1} = body;
            end
            contactList = {};

            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    for k = 1: this.collider.collisions{j}.contactNum
                        if(this.collider.collisions{j}.ground)
                            contact.body1 = this.collider.collisions{j}.constraints{k}.body.index;
                            contact.body2 = 0;
                            contact.normal1 = se3.qRotInv(this.collider.collisions{j}.constraints{k}.body.x0(1:4), this.collider.collisions{j}.constraints{k}.nw);
                            contact.normal2 = [];
                            contact.r1 = this.collider.collisions{j}.constraints{k}.xl;
                            contact.r2 = [];
                        else
                            contact.body1 = this.collider.collisions{j}.constraints{k}.body1.index;
                            contact.body2 = this.collider.collisions{j}.constraints{k}.body2.index;
                            contact.normal1 = se3.qRotInv(this.collider.collisions{j}.constraints{k}.body1.x0(1:4), this.collider.collisions{j}.constraints{k}.nw);
                            contact.normal2 = se3.qRotInv(this.collider.collisions{j}.constraints{k}.body2.x0(1:4), this.collider.collisions{j}.constraints{k}.nw);
                            contact.r1 = this.collider.collisions{j}.constraints{k}.x1;
                            contact.r2 = this.collider.collisions{j}.constraints{k}.x2;
                        end
                        contactList{end+1}=contact;
                    end
                end
            end

            collisionList = {};
            for i = 1 : length(this.collider.activeCollisions)
                for j = this.collider.activeCollisions{i}
                    collision.body1 = this.collider.collisions{j}.body1.index;
                    collision.body2 = this.collider.collisions{j}.body2.index;
                    collision.J1I = this.collider.collisions{j}.J1I;
                    collision.J2I = this.collider.collisions{j}.J2I;
                    collision.mIndices = this.collider.collisions{j}.mIndces;
                    collisionList{end+1} = collision;
                end
            end

            %{
            n = length(contactList);
            Afill = zeros(n,n);
            for i = 1:length(contactList)
                    body1 = bodyList{contactList{i}.body1};
                    I1(1:3) = body1.I;
                    I1(4:6) = body1.m;
                    I1inv = diag(1./I1);
               
                    if(contactList{i}.body2 ~= 0)
                        body2 = bodyList{contactList{i}.body2};
                        I2(1:3) = body2.I;
                        I2(4:6) = body2.m; 
                        I2inv = diag(1./I2);
                    end
            
                for j = i:length(contactList)
                    if(contactList{i}.body1 == contactList{j}.body1)
                        Afill(i,j) = Afill(i,j) + contactList{i}.normal1'*[se3.brac(contactList{i}.r1)', eye(3)] * I1inv *[se3.brac(contactList{j}.r1)', eye(3)]'*contactList{j}.normal1;
                    end
                    if(contactList{i}.body1 == contactList{j}.body2)
                        Afill(i,j) = Afill(i,j) - contactList{i}.normal1'*[se3.brac(contactList{i}.r1)', eye(3)] * I1inv *[se3.brac(contactList{j}.r2)', eye(3)]'*contactList{j}.normal2;
                    end
            
                    if(contactList{i}.body2 ~=0 && contactList{i}.body2 == contactList{j}.body1)
                        Afill(i,j) = Afill(i,j) - contactList{i}.normal2'*[se3.brac(contactList{i}.r2)', eye(3)] * I2inv *[se3.brac(contactList{j}.r1)', eye(3)]'*contactList{j}.normal1;
                    end
                    if(contactList{i}.body2 ~=0 && contactList{i}.body2 == contactList{j}.body2)
                        Afill(i,j) = Afill(i,j) + contactList{i}.normal2'*[se3.brac(contactList{i}.r2)', eye(3)] * I2inv *[se3.brac(contactList{j}.r2)', eye(3)]'*contactList{j}.normal2;
                    end
                end
            end
            Afill = Afill + triu(Afill,1)';
            %}
            file_name = split(this.name);
            save(strcat(file_name{end}, '_states' ,'.mat'),"A","b","contactList","bodyList", "collisionList");
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
