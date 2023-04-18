
classdef (Abstract) Body < handle
	%Body A rigid body connected through a joint to a parent
	
	%%
	properties
		name     % Optional name
		color    % Color for rendering
		density  % Mass/volume
		damping  % Viscous damping
		collide  % Whether collision is enabled
		I_i      % Inertia at body center
        I_inv      % Inertia at body center
        mass     % Mass
        x        % World Position
        x_prev   % Previous World Position
        v        % Velocity
        q        % Quaternion
        q_prev   % Previous quaternion
        w        % Angular Velocity
        E_wi     % Where the body is wrt world
        next     % Next body in traversal order
    end

    	%%
	methods
		%%
		function this = Body(density)
			this.name = ['body'];
			this.color = [0.8 0.8 0.3];
			this.density = density;
			this.damping = 0;
			this.collide = false;
            this.x = [0, 0, 0]';
            this.v = [0, 0, 0]';
            this.q = quaternion(1, 0, 0, 0);
            this.w = [0, 0, 0]';
            this.E_wi = eye(4);

		end
		
		%%
		function setBodyTransform(this, x, q)
			% Sets the transform of this body wrt parent joint
			this.x = x;
			this.q = q;
        end
		
		%%
		function computeInertia(this)
			% Computes inertia at body and joint
			this.computeInertia_();
        end
		%%
		function [Fs,Vs] = draw(this,Fs,Vs)
			if nargin == 1
				Fs = {};
				Vs = {};
			end
			%s = this.getAxisSize();
			%if s > 0
			%	se3.drawAxis(this.E_wi,s);
			%end
			%text(this.E_wi(1,4),this.E_wi(2,4),this.E_wi(3,4),this.name);
			[Fs{end+1},Vs{end+1}] = this.draw_();
			% Go to the next body
			if ~isempty(this.next)
				[Fs,Vs] = this.next.draw(Fs,Vs);
			end
		end
		
		%%
		function s = getAxisSize(this) %#ok<MANU>
			s = 1;
        end
        %%
		function collisions = collideGround(this,groundE,planar,collisions)
			% Runs collision detection with the ground
			if this.collide
				collisions = this.collideGround_(groundE,planar,collisions);
			end
			% Go to the next body
			if ~isempty(this.next)
				collisions = this.next.collideGround(groundE,planar,collisions);
			end
		end

		%%
		function collisions = collideBody(this,that,planar,collisions)
			if this.collide && that.collide
				collisions = this.collideBody_(that,planar,collisions);
			end
		end
        %%
        function updateWithoutConstraints(this, force_ext, h)
            this.x_prev = this.x;
            this.v = this.v + h * force_ext / this.mass;
            this.x = this.x + h * this.v;

            this.q_prev = this.q;
            this.q = this.q + h * 0.5 * quaternion([0 this.w']) * this.q;
            this.q = normalize(this.q);
        end
        %%

        function updateAfterSolve(this, h)
            this.v = (this.x - this.x_prev) / h;
            qinv = quatinv(this.q_prev);
            delta_q = this.q * qinv;
            [qA,qB,qC,qD] = parts(delta_q);
            this.w = 2 * [qB qC qD]' / h;
            if qA < 0
                this.w = -this.w;
            end
            this.E_wi(1:3,1:3) = quat2rotm(this.q);
            this.E_wi(1:3,4) = this.x;
        end

    end
end