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
                    this.solveConSP();
                    this.solveConGS();
					this.t = this.t + this.hs;
					this.ks = this.ks + 1;
				end
				this.k = this.k + 1;
				this.draw();
			end
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
        end

		%%
		function solveConGS(this)

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
