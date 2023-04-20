classdef Scene < handle
	%Scene Test scenes for redmax
	
	properties
		name % scene name
		bodies % list of bodies
		joints % list of joints
		forces % list of forces
		equalities % list of equality constraints
		collisions % list of collisions
		limits % list of joint limits
        constraints %list of constraints
		ground % ground transform, with Z up
		tEnd % end time
		qInit % initial positions
		qdotInit % initial velocities
		t % current time
		h % time step
		k % current step
		mu % coefficient of friction
		restitution % coefficient of restitution
		baumgarte % Baumgarte stabilization
		compliance % compliance
		regDQP % DQP regularization
		planar % 2D planar scene
		Tinit % initial kinetic energy
		Vinit % initial potential energy
		history % step history: does not include the initial state
		auxiliary % auxiliary state for the current step (M, J, etc.)
		nsteps % number of steps to take
        sub_steps% number of substeps per step
		grav % gravity
		drawHz % refresh rate (0 for no draw)
		view % initial viewing angle
		waxis % initial axis
		video %
		computeH % whether to compute the energy
		plotH % whether to plot the energy at the end
		Hexpected % expected energy
		task % task for the adjoint method
		solverInfo % solver info for the current step
	end
	
	methods
		function this = Scene()
			% Default values
			this.name = '';
			this.bodies = {};
			this.joints = {};
			this.equalities = {};
			this.collisions = {};
			this.limits = {};
            this.constraints = {};
			this.ground.E = zeros(4);
			this.ground.size = 10;
			this.tEnd = 1;
			this.qInit = [];
			this.qdotInit = [];
			this.h = 1e-2;
			this.t = 0;
			this.k = 0;
			this.mu = 0;
			this.restitution = 0;
			this.baumgarte = 0;
			this.compliance = 0;
			this.regDQP = 1e-4;
			this.planar = false;
			this.Tinit = 0;
			this.Vinit = 0;
			this.history = [];
			this.auxiliary = [];
			this.nsteps = 0;
			this.grav = [0 0 -980]';
			this.drawHz = 15;
			this.view = 3;
			this.waxis = [];
			this.video = [];
			this.computeH = true;
			this.plotH = true;
			this.Hexpected = zeros(1,2);
			this.task = [];
			this.solverInfo = [];
		end
		
		%%
		function init(this)
			if usejava('jvm')
				colormap('default'); % Restore colormap mode
            end

			% Initialize bodies
			nbodies = length(this.bodies);
			for i = 1 : nbodies
				this.bodies{i}.computeInertia();
				if i < nbodies
					this.bodies{i}.next = this.bodies{i+1}; %#ok<*SAGROW>
				end
			end
						
			% Other initial values
			this.nsteps = ceil(this.tEnd/this.h);
            this.sub_steps = 30;
        end

        %%
        function solvePositions(this)
            for i = 1 : length(this.constraints)
                this.constraints{i}.solvePositions();
            end
            %{
            for i = 1 : length(this.collisions)
				collision = this.collisions{i};
				nbodies = length(collision.bodies);
                % Ground collisions 
				if nbodies == 1
                    body = collision.bodies{1};
                    r1 = body.E_wi * [collision.xl{1}; 1];
                    r1 = r1(1:3);
                    n = collision.nw;
                    c = collision.d;
                    temp = cross(r1, n);
                    w1 = 1 / body.mass + temp' * body.I_inv * temp;
                    delta_lambda = -c / w1;
                    p = delta_lambda * n;
                    body.x = body.x + p / body.mass;
                    temp = body.I_inv * cross(r1, p);
                    body.q = body.q + 0.5 * quaternion([0 temp']) * body.q;
                end
            end
            %}
		end

		%%
		function collide(this)
			% Clear collisions
			this.collisions = {};
			
			% Ground collisions (only if ground.E isn't the zero matrix)
			if this.ground.E(4,4) == 1
				this.collisions = this.bodies{1}.collideGround(this.ground.E,this.planar,this.collisions);
			end

			% Body-body collisions
			for i = 1 : length(this.bodies)
				bodyI = this.bodies{i};
				for j = i+1 : length(this.bodies)
					bodyJ = this.bodies{j};
					% Are bodyI and bodyJ parent/child?
					if ~isempty(bodyI.joint.parent)
						if bodyI.joint.parent.body == bodyJ
							continue;
						end
					end
					if ~isempty(bodyJ.joint.parent)
						if bodyJ.joint.parent.body == bodyI
							continue;
						end
					end
					% Check for collisions
					this.collisions = bodyI.collideBody(bodyJ,this.planar,this.collisions);
				end
			end
			% Look for collisions between a free body and a fixed body.
			% These will be treated as ground collisions.
			indicesToKeep = [];
			for i = 1 : length(this.collisions)
				collision = this.collisions{i};
				nbodies = length(collision.bodies);
				if nbodies == 1
					ndof = 6;
					if ndof == 0
						% Completely remove this constraint
					else
						% This body is free, so keep this constraint
						indicesToKeep(end+1) = i; %#ok<AGROW> 
					end
				else
					ndof0 = length(collision.bodies{1}.joint.getAncestorIndicesR());
					ndof1 = length(collision.bodies{2}.joint.getAncestorIndicesR());
					if ndof0 == 0 && ndof1 == 0
						% Completely remove this constraint
					elseif ndof0 == 0 && ndof1 > 0
						% Remove first body from this constraint
						indicesToKeep(end+1) = i; %#ok<AGROW> 
						this.collisions{i}.bodies = {this.collisions{i}.bodies{2}}; %#ok<CCAT1> 
						this.collisions{i}.xl = {this.collisions{i}.xl{2}}; %#ok<CCAT1> 
						% The normal should be wrt the new first body
						this.collisions{i}.nw = -this.collisions{i}.nw;
					elseif ndof0 > 0 && ndof1 == 0
						% Remove second body from this constraint
						indicesToKeep(end+1) = i; %#ok<AGROW> 
						this.collisions{i}.bodies = {this.collisions{i}.bodies{1}}; %#ok<CCAT1> 
						this.collisions{i}.xl = {this.collisions{i}.xl{1}}; %#ok<CCAT1> 
					else
						% Both bodies are free, so keep this constraint
						indicesToKeep(end+1) = i; %#ok<AGROW> 
					end
				end
			end
			this.collisions = this.collisions(indicesToKeep);
        end

        %%
		function draw(this)
			if this.t == 0
				clf;
				xlabel('X');
				ylabel('Y');
				zlabel('Z');
				axis equal;
				if ~isempty(this.waxis)
					axis(this.waxis);
				end
				ax = gca;
				ax.Clipping = 'off';
				grid on;
				view(this.view); %#ok<CPROP>
			end
			if this.drawHz > 0 && (floor(this.t*this.drawHz) > floor((this.t-this.h)*this.drawHz))
				cla;
				hold on;

				[Fs,Vs] = this.bodies{1}.draw();
				if ~isempty(this.task)
					this.task.draw();
				end

				% Draw ground
				se3.drawAxis(this.ground.E);
				s = this.ground.size/2;
				V = this.ground.E(1:3,:)*[-s -s 0 1; s -s 0 1; s s 0 1; -s s 0 1]';
				F = [1 2 3 4];
				patch('Faces',F,'Vertices',V','FaceColor',[0.9 0.9 0.9]);

				% Lighting
				l = light('Style','local','Position',[0 -50 100]);

				% Shadow:
				%   p = l - ((d + dot(n,l))/dot(n,v - l))*(v - l),
				% where l is the light position, {n,d} is the plane, and v is the
				% vertex position.
				if this.ground.E(4,4) ~= 0
					n = this.ground.E(1:3,3);
					d = -n'*this.ground.E(1:3,4);
					lpos = l.Position';
					for i = 1 : length(Fs)
						F = Fs{i};
						V = Vs{i}';
						Vl = V - lpos;
						Vshadow = lpos - ((d + n'*lpos)./(n'*Vl)).*Vl + 1e-5*n;
						patch('Faces',F,'Vertices',Vshadow','EdgeColor','none','FaceColor',[0.2 0.2 0.2]);
					end
				end

				title(sprintf('t = %.4f',this.t));
				drawnow;

				if ~isempty(this.video)
					this.video.writeVideo(getframe(gcf));
				end
			end
		end
    end
end

