classdef BodyCuboid < Body
	%%
	properties
		sides       % Side lengths
		collision2D % Apply 2D collision along the X=0, Y=0, or Z=0 plane
	end
	
	%%
	methods
		%%
		function this = BodyCuboid(density,sides)
			this = this@Body(density);
			this.sides = sides;
			this.collision2D = '';
		end
		
		%%
		function computeInertia_(this)
			% Computes inertia at body
            temp = se3.inertiaCuboid(this.sides,this.density);
			this.I_i = diag(temp(1:3));
            this.mass = temp(4);
            this.I_inv = inv(this.I_i);
		end

		%%
		function collisions = collideGround_(this,groundE,planar,collisions)
			% Collision detection
			xg = groundE(1:3,4); % ground origin
			ng = groundE(1:3,3); % ground normal
			S = eye(4);
			S(1:3,1:3) = diag(0.5*this.sides);
			switch this.collision2D
				case 'X'
					xl = S*[
						 0 -1 -1  1
						 0 -1  1  1
						 0  1 -1  1
						 0  1  1  1
						]';
				case 'Y'
					xl = S*[
						-1  0 -1  1
						-1  0  1  1
						 1  0 -1  1
						 1  0  1  1
						]';
				case 'Z'
					xl = S*[
						-1 -1  0  1
						-1  1  0  1
						 1 -1  0  1
						 1  1  0  1
						]';
				otherwise
					xl = S*[
						-1 -1 -1  1
						-1 -1  1  1
						-1  1 -1  1
						-1  1  1  1
						 1 -1 -1  1
						 1 -1  1  1
						 1  1 -1  1
						 1  1  1  1
						]';
			end
			%xl = xl(:,1); % for debugging JointPlanar
			%xl = xl(:,end); % for debugging
			%xl = xl(:,[1,5]); % for debugging
			xw = this.E_wi*xl;
			for i = 1 : size(xl,2)
				xli = xl(1:3,i);
				xwi = xw(1:3,i);
				% penetration depth
				d = ng'*(xwi - xg);
				if d > 0
					% No collision
					continue
				end
				
				% Contact point world velocity
				collision.bodies{1} = this;
				collision.xl{1} = xli; % local pos
				collision.xw = xwi; % world pos
				collision.nw = ng; % world nor from ground
				%collision.Tw = this.getTangentBasis(collision,planar); % world tangent based on velocity
				collision.d = d; % penetration depth
				%collision.l = []; % initial Lagrange multiplier
				%collision.idxB = []; % index into the global 
				collisions{end+1} = collision; %#ok<AGROW> 
			end
		end

		%%
		function collisions = collideBody_(this,that,planar,collisions)
			if ~isa(that,'redmax.BodyCuboid')
				% Only Cuboid-Cuboid is supported
				return;
			end
			c = odeBoxBox_mex(this.E_wi,this.sides,that.E_wi,that.sides);
			for i = 1 : c.count
				xw = c.pos(:,i);
				if this.collision2D == 'Z'
					% Assume that the 2D scene is on the z=0 plane
					if xw(3) < 0
						% Remove collisions with negative z
						continue;
					end
					% Make the positive collisions happen at z=0
					xw(3) = 0;
				end
				nw = -c.nor; % sign depends on the convention of the collision detector
				d = -c.depth(i);
				collision.bodies{1} = this;
				collision.bodies{2} = that;
				collision.xl{1} = this.E_iw(1:3,:)*[xw;1];
				collision.xl{2} = that.E_iw(1:3,:)*[xw;1];
				collision.xw = xw;
				collision.nw = nw;
				collision.Tw = this.getTangentBasis(collision,planar); % world tangent based on velocity
				collision.d = d;
				collisions{end+1} = collision; %#ok<AGROW> 
			end
		end
		
		%%
		function [F,V] = draw_(this)
			[F,V] = se3.patchCuboid(this.E_wi,this.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',this.color);
		end
		
		%%
		function s = getAxisSize(this)
			s = min(this.sides);
		end
	end
end