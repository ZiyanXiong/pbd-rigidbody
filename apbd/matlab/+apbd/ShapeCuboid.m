classdef ShapeCuboid < apbd.Shape
	%%
	properties
		sides % Side lengths
	end
	
	%%
	methods
		%%
		function this = ShapeCuboid(sides)
			this = this@apbd.Shape();
			this.sides = sides;
		end
		
		%%
		function I = computeInertia(this,density)
            % I is the 6x1 diagonal rigid inertia, assuming that the frame origin is at the
            % center of mass and the axes are oriented along the principal axes. We store
            % the rotation on top of translations, so that I(1:3) is the rotational inertia
            % and I(4:6) is the translational inertia.
            I = se3.inertiaCuboid(this.sides,density);
		end

		%%
		function flag = broadphaseGround(this,E,Eg)
			% Check the height of the center
			xl = [0 0 0 1]';
			xw = E*xl;
			xg = Eg\xw;
			r = norm(this.sides/2); % dist to a corner
			flag = xg(3) < 1.5*r;
		end

		%%
		function cdata = narrowphaseGround(this,E,Eg)
			cdata = [];
			s = this.sides/2;
			xl = ones(4,8); % local space
			xl(1:3,1) = [-s(1) + 0.0, -s(2), -s(3)]';
			xl(1:3,2) = [ s(1), -s(2), -s(3)]';
			xl(1:3,3) = [ s(1),  s(2), -s(3)]';
			xl(1:3,4) = [-s(1) + 0.0,  s(2), -s(3)]';
			xl(1:3,5) = [-s(1) + 0.0, -s(2),  s(3)]';
			xl(1:3,6) = [ s(1), -s(2),  s(3)]';
			xl(1:3,7) = [ s(1),  s(2),  s(3)]';
			xl(1:3,8) = [-s(1) + 0.0,  s(2),  s(3)]';
			% Debug
			%xl(1,[1,4,5,8]) = -s(1)*0.5;
			%xl(1,[2,3,6,7]) = s(1)*0.5;
			% End debug
			xw = E*xl;
			xg = Eg\xw;
			for i = 1 : 8
				% This only supports vertex collisions
				d = xg(3,i);
				if d < 0
					xgproj = xg(:,i);
					xgproj(3) = 0; % project onto the floor plane
					cdata(end+1).d = d; %#ok<AGROW>
					cdata(end).xl = xl(1:3,i);
					cdata(end).xw = Eg(1:3,:)*xgproj; % transform to world space
					cdata(end).nw = Eg(1:3,3); % normal
					cdata(end).vw = [0 0 0]'; % unused for now, assming the ground doesn't move
				end
			end
		end

		%%
		function flag = broadphaseShape(this,E1,that,E2)
			if isa(that,'apbd.ShapeCuboid')
				p1 = E1(1:3,4);
				p2 = E2(1:3,4);
				d = norm(p1 - p2);
				r1 = norm(this.sides/2); % dist to a corner
				r2 = norm(that.sides/2); % dist to a corner
				flag = d < 1.5*(r1 + r2);
			else
				error('Unsupported shape');
			end
		end

		%%
		function cdata = narrowphaseShape(this,E1,that,E2)
			cdata = [];
			if isa(that,'apbd.ShapeCuboid')
				R1 = E1(1:3,1:3);
				R2 = E2(1:3,1:3);
				p1 = E1(1:3,4);
				p2 = E2(1:3,4);
				s1 = this.sides;
				s2 = that.sides;
				collisions = odeBoxBox_mex(E1,s1,E2,s2);
				nw = collisions.nor; % The normal is outward from body 1 (red)
				n1 =  R1'*nw;
				n2 = -R2'*nw; % negate since nw is defined wrt body 1
				for i = 1 : collisions.count
					xw = collisions.pos(:,i);
					d = -collisions.depth(i); % odeBoxBox returns positive depth for hits
					% Compute local point on body 1 with ray casting
					x1 = R1'*(xw - p1);
					[~,t1] = this.raycast(x1,n1);
					x1 = x1 - t1*n1; % negate since smits_mul returns negative t for rays starting inside the box
					% Compute local point on body 2 with ray casting
					x2 = R2'*(xw - p2);
					[~,t2] = that.raycast(x2,n2);
					x2 = x2 - t2*n2; % negate since smits_mul returns negative t for rays starting inside the box
					cdata(end+1).d = d; %#ok<AGROW>
					cdata(end).xw = xw;
					cdata(end).nw = nw;
					cdata(end).x1 = x1;
					cdata(end).x2 = x2;
				end
			else
				error('Unsupported shape');
			end
		end

		%%
		function [hit,t] = raycast(this,x,n)
			thresh = 1e-6;
			bmax = 0.5*this.sides';
			bmin = -bmax;
			x = (1 - thresh)*x; % make the point go slightly inside the box
			[hit,t] = apbd.ShapeCuboid.smits_mul(x,-n,bmin,bmax); % negate ray since it starts inside the box
		end


		%%
		function [F,V] = draw(this,E,color,axisSize)
			if nargin < 4
				axisSize = 0;
			end
			if nargin < 3
				color = [0.5 0.5 0.5];
			end
			[F,V] = se3.patchCuboid(E,this.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',color);
			if axisSize > 0
				se3.drawAxis(E,axisSize);
			end
		end
	end

	methods (Static)
		% /******************************************************************************
		%
		%   This source code accompanies the Journal of Graphics Tools paper:
		%
		%   "Fast Ray-Axis Aligned Bounding Box Overlap Tests With Pluecker Coordinates" by
		%   Jeffrey Mahovsky and Brian Wyvill
		%   Department of Computer Science, University of Calgary
		%
		%   This source code is public domain, but please mention us if you use it.
		%
		%  ******************************************************************************/
		function [hit,t] = smits_mul(rx, rv, bmin, bmax)
			% rx: ray origin (3x1)
			% rv: ray direction (3x1)
			% bmin: box min (3x1)
			% bmax: box max (3x1)

			r = apbd.ShapeCuboid.make_ray(rx(1), rx(2), rx(3), rv(1), rv(2), rv(3));
			b = apbd.ShapeCuboid.make_aabox(bmin(1), bmin(2), bmin(3), bmax(1), bmax(2), bmax(3));

			tnear = -1e6;
			tfar = 1e6;

			hit = false;
			t = 0;

			% multiply by the inverse instead of dividing
			t1 = (b.x0 - r.x) * r.ii;
			t2 = (b.x1 - r.x) * r.ii;

			if t1 > t2
				temp = t1;
				t1 = t2;
				t2 = temp;
			end
			if t1 > tnear
				tnear = t1;
			end
			if t2 < tfar
				tfar = t2;
			end
			if tnear > tfar
				return;
			end
			if tfar < 0.0
				return;
			end

			t1 = (b.y0 - r.y) * r.ij;
			t2 = (b.y1 - r.y) * r.ij;

			if t1 > t2
				temp = t1;
				t1 = t2;
				t2 = temp;
			end
			if t1 > tnear
				tnear = t1;
			end
			if t2 < tfar
				tfar = t2;
			end
			if tnear > tfar
				return;
			end
			if tfar < 0.0
				return
			end

			t1 = (b.z0 - r.z) * r.ik;
			t2 = (b.z1 - r.z) * r.ik;

			if t1 > t2
				temp = t1;
				t1 = t2;
				t2 = temp;
			end
			if t1 > tnear
				tnear = t1;
			end
			if t2 < tfar
				tfar = t2;
			end
			if tnear > tfar
				return
			end
			if tfar < 0.0
				return
			end

			t = tnear;
			hit = true;

		end

		%%
		function r = make_ray(x, y, z, i, j, k)
			r.x = x;
			r.y = y;
			r.z = z;
			r.i = i;
			r.j = j;
			r.k = k;
			r.ii = 1.0 / i;
			r.ij = 1.0 / j;
			r.ik = 1.0 / k;
		end

		%%
		function a = make_aabox(x0, y0, z0, x1, y1, z1)
			if x0 > x1
				a.x0 = x1;
				a.x1 = x0;
			else
				a.x0 = x0;
				a.x1 = x1;
			end
			if y0 > y1
				a.y0 = y1;
				a.y1 = y0;
			else
				a.y0 = y0;
				a.y1 = y1;
			end
			if z0 > z1
				a.z0 = z1;
				a.z1 = z0;
			else
				a.z0 = z0;
				a.z1 = z1;
			end
		end

	end

end
