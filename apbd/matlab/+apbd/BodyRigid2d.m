classdef BodyRigid2d < apbd.Body
	%BodyRigid2d A 2D rigid body in the XY plane
	% This class is not faster than 3D but is useful for debugging.
	
	%%
	properties
		shape    % Associated shape
		density  % Mass/volume
		Mr       % Rotational inertia (3x1)
		Mp       % Translational inertia (1x1)

		% Drawing etc.
		color    % Color for rendering
		axisSize % 
    end

	%%
	methods
		%%
		function this = BodyRigid2d(shape,density)
			global CM; %#ok<GVMIS>

			this = this@apbd.Body(4);
			this.shape = shape;
			this.density = density;
			this.Mr = zeros(3,1);
			this.Mp = 0;

			this.color = CM(mod(this.index-1,size(CM,1))+1,:);
			this.axisSize = 1;
		end
		
		%%
		function setInitTransform(this,E)
			q = se3.matToQ(E(1:3,1:3));
			this.xInit(1:2) = q(3:4); % just the zw component
			this.xInit(3:4) = E(1:2,4); % just the xy components
		end

		%%
		function Einit = computeInitTransform(this)
			% This can be used only after calling setInitTransform()
			Einit = eye(4);
			[q,p] = apbd.BodyRigid2d.unproj(this.xInit);
			Einit(1:3,1:3) = se3.qToMat(q);
            Einit(1:2,4) = p;
		end

		%%
		function setInitVelocity(this,phi)
			q = apbd.BodyRigid2d.unproj(this.xInit);
			pdot = se3.qRot(q,phi(4:6));
			this.xdotInit(3:4) = pdot(1:2);
			qdot = se3.wToQdot(q,phi(1:3));
			this.xdotInit(1:2) = qdot(3:4);
		end
		
		%%
		function E = computeTransform(this)
			E = eye(4);
			[q,p] = apbd.BodyRigid2d.unproj(this.x);
			E(1:3,1:3) = se3.qToMat(q);
            E(1:3,4) = p;
		end

		%%
		function xw = transformPoint(this,xl)
			[q,p] = apbd.BodyRigid2d.unproj(this.x);
			xw = se3.qRot(q,xl) + p;
		end

		%%
		function vw = transformVector(this,vl)
			q = apbd.BodyRigid2d.unproj(this.x);
			vw = se3.qRot(q,vl);
		end

		%%
		function xl = invTransformPoint(this,xw)
			[q,p] = apbd.BodyRigid2d.unproj(this.x);
			xl = se3.qRotInv(q,xw - p);
		end

		%%
		function vl = invTransformVector(this,vw)
			q = apbd.BodyRigid2d.unproj(this.x);
			vl = se3.qRotInv(q,vw);
		end
		
		%%
		function computeInertiaConst(this)
			% Computes inertia for the shape
			I = this.shape.computeInertia(this.density);
			this.Mr = I(1:3);
            this.Mp = I(4);
		end

		%%
		function v = computePointVel(this,xl,k,ks,hs)
			xdot = this.computeVelocity(k,ks,hs);
			[qdot,pdot] = apbd.BodyRigid2d.unprojVel(xdot);
			q = apbd.BodyRigid2d.unproj(this.x);
			w = se3.qdotToW(q,qdot); % angular velocity in body coords
			% v = R*cross(w,xl) + pdot
			v = se3.qRot(q,se3.cross(w,xl)) + pdot;
		end
		
        %%
		function stepBDF1(this,k,ks,hs,grav)
			xdot = this.computeVelocity(k,ks,hs);
			[qdot,v] = apbd.BodyRigid2d.unprojVel(xdot);
            this.x0 = this.x;
			[q,p] = apbd.BodyRigid2d.unproj(this.x);
			w = se3.qdotToW(q,qdot); % angular velocity in body coords
			f = zeros(3,1); % translational force in world space
			t = zeros(3,1); % angular torque in body space
			m = this.Mp; % scalar mass
			I = this.Mr; % inertia in body space
			Iw = I.*w; % angular momentum in body space
			f = f + m*grav; % Gravity
			t = t + se3.cross(Iw,w); % Coriolis
			% Integrate velocities
			w = w + hs*(I.\t);
			v = v + hs*(m \f);
			qdot = se3.wToQdot(q,w);
			% Integrate positions
			q = q + hs*qdot;
			p = p + hs*v;
			q = q/norm(q);
			%fprintf('%f %f\n',q(3:4));
			this.x(1:2) = q(3:4);
			this.x(3:4) = p(1:2);
		end

		%%
		function [T,V] = computeEnergies(this,k,ks,hs,grav)
			xdot = this.computeVelocity(k,ks,hs);
			[qdot,v] = apbd.BodyRigid2d.unprojVel(xdot);
			[q,p] = apbd.BodyRigid2d.unproj(this.x);
			w = se3.qdotToW(q,qdot); % angular velocity in body coords
			m = this.Mp; % scalar mass
			I = this.Mr; % inertia in body space
			Iw = I.*w; % angular momentum in body space
			% Energy
			T = 0.5*w'*Iw + 0.5*m*(v'*v);
			V = -m*grav'*p;
		end

		%%
		function flag = broadphaseGround(this,Eg)
			E = this.computeTransform();
			flag = this.shape.broadphaseGround(E,Eg);
		end

		%%
		function cdata = narrowphaseGround(this,Eg)
			E = this.computeTransform();
			cdata = this.shape.narrowphaseGround(E,Eg);
			% Remove all duplicate collisions along the z direction
			keep = [];
			for k = 1 : length(cdata)
				c = cdata(k);
				if c.xl(3) < 0
					keep(end+1) = k; %#ok<AGROW>
					cdata(k).xl(3) = 0;
					cdata(k).xw(3) = 0;
				end
			end
			cdata = cdata(keep);
		end

		%%
		function flag = broadphaseRigid(this,that)
			E1 = this.computeTransform();
			E2 = that.computeTransform();
			flag = this.shape.broadphaseShape(E1,that.shape,E2); % bad syntax
		end

		%%
		function cdata = narrowphaseRigid(this,that)
			E1 = this.computeTransform();
			E2 = that.computeTransform();
			cdata = this.shape.narrowphaseShape(E1,that.shape,E2); % bad syntax
			% Remove all duplicate collisions along the z direction
			keep = [];
			for k = 1 : length(cdata)
				c = cdata(k);
				if c.xl(3) < 0
					keep(end+1) = k; %#ok<AGROW>
					cdata(k).xl(3) = 0;
					cdata(k).xw(3) = 0;
				end
			end
			cdata = cdata(keep);
		end

		%%
		function [F,V] = draw(this)
			E = this.computeTransform();
			[F,V] = this.shape.draw(E,this.color,this.axisSize);
		end
	end

	methods (Static)
		%%
		function [q,p] = unproj(x)
			q = zeros(4,1);
			q(3:4) = x(1:2);
			if nargout == 2
				p = zeros(3,1);
				p(1:2) = x(3:4);
			end
		end

		%%
		function [qdot,pdot] = unprojVel(xdot)
			qdot = zeros(4,1);
			qdot(3:4) = xdot(1:2);
			if nargout == 2
				pdot = zeros(3,1);
				pdot(1:2) = xdot(3:4);
			end
		end

		%%
		function J = jacobian(xl,w)
			wx = w(1);
			wy = w(2);
			wz = w(3);
			w0 = w(4);
			Jaw = 2*[
				 wx, -wy, -wz,  w0
				 wy,  wx, -w0, -wz
				 wz,  w0,  wx,  wy
				 wy,  wx,  w0,  wz
				-wx,  wy, -wz,  w0
				-w0,  wz,  wy, -wx
				 wz, -w0,  wx, -wy
				 w0,  wz,  wy,  wx
				-wx, -wy,  wz,  w0
				];
			Jxa = apbd.BodyAffine.jacobian(xl);
			J = [Jxa(:,1:9)*Jaw, eye(3)];
		end
	end
end
