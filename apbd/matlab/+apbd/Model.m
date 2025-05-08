classdef Model < handle
	%Model Class to hold model data

	properties
		bodies % list of bodies
		constraints %list of constraints
        joints %list of joints
		collider % collision handler
		grav % gravity
        f % force for each body at current timestep
        tau % torque for each body at current timestep
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
        solverType % 1: TGS 2: 2PSP
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
            this.joints = {};
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
			this.k = 1;
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
			while this.k <= this.steps
				this.ks = 0;
                if this.k == 20
                    fprintf("Pause.");
                end
                %{
                if(this.solverType == 1)
                    this.solveConTGS();
                else
                    this.solveConGlobal();
                end
                %}
                this.solveConGlobal();
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
            this.f = zeros(length(this.bodies), 3);
            this.tau = zeros(length(this.bodies), 3);
            for i = 1 : length(this.joints)
                [this.f, this.tau] = this.joints{i}.applyForceTorque(this.f, this.tau, this.k);
            end
			for i = 1 : length(this.bodies)
			    this.bodies{i}.stepBDF1(this.h,this.grav,this.f(i,:)',this.tau(i,:)');
			end
        end

		%%
		function solveConTGS(this)
            this.collider.run();
            for i = this.collider.activeCollisions
                this.collider.collisions{i}.initConstraints(this.h, this.hs);
            end
            for i = 1 : length(this.joints)
                this.joints{i}.init(this.h, this.hs);
            end
            this.draw();

			this.stepBDF1();
            for i = this.collider.activeCollisions
                this.collider.collisions{i}.compute_b0();
            end

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
                for i = this.collider.activeCollisions
                    this.collider.collisions{i}.solveCollisionNor(false);
                    this.collider.collisions{i}.solveCollisionTan(false);
                end

                % Gauss-Seidal 
                for i = 1 : length(this.joints)
                    this.joints{i}.solveCollisionNor(false);
                    this.joints{i}.solveCollisionTan(false);
                end

                for i = 1 : length(this.bodies)
                    this.bodies{i}.updateStates(this.hs);
                end
                %this.draw();
			    this.t = this.t + this.hs;
			    this.ks = this.ks + 1;
            end
            this.iterVec(end+1) = this.substeps;

            for i = 1 : length(this.bodies)
                this.bodies{i}.initVelocitySolve(this.h);
            end

            for iter = 1: 10
                for i = this.collider.activeCollisions
                    this.collider.collisions{i}.solveCollisionNorVel(10);
                    this.collider.collisions{i}.solveCollisionTanVel(10);
                end
                for i = 1 : length(this.bodies)
                    this.bodies{i}.updateVelocities(1/10);
                end
            end

            for i = 1 : length(this.bodies)
                this.bodies{i}.integrateStates();
            end

            n = 1;
            ci = 1;
            for i = this.collider.activeCollisions
                this.collider.collisions{i}.index = ci;
                this.collider.collisions{i}.mIndces = n : n - 1 + this.collider.collisions{i}.contactNum * 3;
                n = n + this.collider.collisions{i}.contactNum * 3;
                this.collider.collisions{i}.compute_b();
                ci = ci + 1;
            end
            rs = zeros(n-1,1);
            for i = this.collider.activeCollisions
                rs(this.collider.collisions{i}.mIndces) = this.collider.collisions{i}.b;
            end

            this.rVec(end+1) = norm(rs(rs>0));
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


		%%
        function solveConGlobal(this)
            this.collider.run();
            for i = this.collider.activeCollisions
                this.collider.collisions{i}.initConstraints(this.h, this.hs);
            end
            for i = 1 : length(this.joints)
                this.joints{i}.init(this.h, this.hs, this.k);
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
                for i = this.collider.activeCollisions
                    this.collider.collisions{i}.index = ci;
                    this.collider.collisions{i}.mIndces = n : n - 1 + this.collider.collisions{i}.contactNum * 3;
                    n = n + this.collider.collisions{i}.contactNum * 3;
                    this.collider.collisions{i}.computeJ_b();
                    this.collider.collisions{i}.compute_d();
                    ci = ci + 1;
                end

                for i = 1 : length(this.joints)
                    this.joints{i}.index = ci;
                    this.joints{i}.mIndces = n : n - 1 + this.joints{i}.lambdaLen;
                    n = n + this.joints{i}.lambdaLen;
                    this.joints{i}.computeJ_b();
                    this.joints{i}.compute_d();
                    ci = ci + 1;
                end

                n = n-1;
                b = zeros(n,1);
                d = zeros(n,1);
                for i = this.collider.activeCollisions
                    inds = this.collider.collisions{i}.mIndces;
                    b(inds) = this.collider.collisions{i}.b;
                    d(inds) = this.collider.collisions{i}.d;
                end

                for i = 1 : length(this.joints)
                    inds = this.joints{i}.mIndces;
                    b(inds) = this.joints{i}.b;
                    d(inds) = this.joints{i}.d;
                end

                for i = 1:length(this.bodies)
                    this.bodies{i}.colIndices = (i - 1)*6+1 : i*6;
                end
                L = zeros(n,6*length(this.bodies));

                for i = this.collider.activeCollisions
                    rows = this.collider.collisions{i}.mIndces;
                    if(this.collider.collisions{i}.ground)
                        cols =  this.collider.collisions{i}.body1.colIndices;
                        L(rows,cols) = this.collider.collisions{i}.J1I;
                    else
                        cols =  this.collider.collisions{i}.body1.colIndices;
                        L(rows,cols) = this.collider.collisions{i}.J1I;
                        cols =  this.collider.collisions{i}.body2.colIndices;
                        L(rows,cols) = this.collider.collisions{i}.J2I;
                    end
                end

                for i = 1 : length(this.joints)
                    rows = this.joints{i}.mIndces;
                    if(this.joints{i}.ground)
                        cols =  this.joints{i}.body1.colIndices;
                        L(rows,cols) = this.joints{i}.J1I;
                    else
                        cols =  this.joints{i}.body1.colIndices;
                        L(rows,cols) = this.joints{i}.J1I;
                        cols =  this.joints{i}.body2.colIndices;
                        L(rows,cols) = this.joints{i}.J2I;
                    end
                end

                A = L*L';
                mu = this.bodies{1}.mu;

                itermax = this.substeps;                
                solver = ConstraintSolver(itermax,1e-6);
                if(isempty(this.joints))
                    contactConstraintEndInd = length(b);
                else
                    contactConstraintEndInd = this.joints{1}.mIndces(1) - 1;
                end
                %[lambdas, lambdav] = solver.Cone_GPQP(A, b, mu, contactConstraintEndInd);
                %[lambdas, lambdav] = solver.Temporal_Gauss_Sidiel(A, b, d, contactConstraintEndInd, mu, 150,lambdas);
                %[lambdas, lambdav] = solver.SOCP(L, b, mu);
                if(this.solverType == 1)
                    [lambdas, lambdav] = solver.Gauss_Sidiel(A, b, mu, contactConstraintEndInd);
                    lambdas = pinv(A)*b;
                else
                    [lambdas, lambdav] = solver.Staggered(A, b, mu, contactConstraintEndInd);
                end
                this.iterVec(end+1) = solver.itercount;
                this.rVec(end+1) = solver.rs(this.iterVec(end));
                if(this.k == 10)
                    if(this.solverType == 1)
                        fileName = sprintf("rVec_GS_step_%d.mat", this.k);
                        rs_gs = solver.rs;
                        save(strcat(this.resultFolder, fileName), 'rs_gs'); 
                    else
                        fileName = sprintf("rVec_GPQP_step_%d.mat", this.k);
                        rs_gpqp = solver.rs;
                        save(strcat(this.resultFolder, fileName), 'rs_gpqp'); 
                    end
                end

                for i = this.collider.activeCollisions
                    l = this.collider.collisions{i}.contactNum*3 - 1;
                    start = this.collider.collisions{i}.mIndces(1);
                    lambdai = lambdas(start:start+l);
                    for j = 1: this.collider.collisions{i}.contactNum
                        this.collider.collisions{i}.constraints{j}.applyLambda(lambdai(3*(j-1) + 1: 3*j));
                    end
                end

                for i = 1 : length(this.joints)
                    l = this.joints{i}.lambdaLen-1;
                    start = this.joints{i}.mIndces(1);
                    lambdai = lambdas(start:start+l);
                    this.joints{i}.applyLambdas(lambdai);
                    this.joints{i}.compute_b();
                    this.joints{i}.recordTorques(this.k);
                end

                for i = 1 : length(this.bodies)
                    this.bodies{i}.updateStatesDirect(this.h);
                end
            end

            %{
            [~,lambdav] = solver.Temporal_Gauss_Sidiel(A, b+d, zeros(n,1), mu, contactConstraintEndInd, 10, lambdav);
            dlambdas = lambdav - lambdas;
            %dlambdavs = lambdavs - lambdas;

            for i = this.collider.activeCollisions
                l = this.collider.collisions{i}.contactNum*3 - 1;
                start = this.collider.collisions{i}.mIndces(1);
                lambdai = dlambdas(start:start+l);
                for j = 1: this.collider.collisions{i}.contactNum
                    this.collider.collisions{i}.constraints{j}.applyLambda(lambdai(3*(j-1) + 1: 3*j));
                end
            end

            for i = 1 : length(this.joints)
                l = this.joints{i}.lambdaLen-1;
                start = this.joints{i}.mIndces(1);
                lambdai = dlambdas(start:start+l);
                this.joints{i}.applyLambdas(lambdai);
                this.joints{i}.compute_b();
            end
            %}
            

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

                for i = 1: length(this.joints)
                    this.joints{i}.draw();
                end

				% Draw collisions
                for i = this.collider.activeCollisions
                    %this.collider.collisions{i}.draw();
                end

				% Lighting
				l = light('Style','local','Position',[50 20 100]);

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
