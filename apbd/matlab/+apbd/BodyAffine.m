classdef BodyAffine < apbd.Body
	%BodyAffine An affine body
	
	%%
	properties
		shape    % Associated shape
		density  % Mass/volume
		Wa       % Inverse affine inertia (3x1)
		Wp       % Inverse translational inertia (1x1)

		% Drawing etc.
		color    % Color for rendering
		axisSize % 
    end

	%%
	methods
		%%
		function this = BodyAffine(shape,density)
			global CM; %#ok<GVMIS>

			this = this@apbd.Body(12);
			this.shape = shape;
			this.density = density;
			this.Wa = zeros(3,1);
			this.Wp = 0;

			this.color = CM(mod(this.index-1,size(CM,1))+1,:);
			this.axisSize = 1;
		end
		
		%%
		function setInitTransform(this,E)
			this.xInit(1:9) = reshape(E(1:3,1:3)',9,1);
			this.xInit(10:12) = E(1:3,4);
		end

		%%
		function Einit = computeInitTransform(this)
			% This can be used only after calling setInitTransform()
			Einit = eye(4);
			Einit(1:3,1:3) = reshape(this.xInit(1:9),3,3)';
            Einit(1:3,4) = this.xInit(10:12);
		end

		%%
		function setInitVelocity(this,phi)
			E = this.computeInitTransform();
			Edot = E*se3.brac(phi);
			this.xdotInit(1:9) = reshape(Edot(1:3,1:3)',9,1);
			this.xdotInit(10:12) = Edot(1:3,4);
		end
		
		%%
		function E = computeTransform(this)
			E = eye(4);
			E(1:3,1:3) = reshape(this.x(1:9),3,3)';
            E(1:3,4) = this.x(10:12);
		end

		%%
		function xw = transformPoint(this,xl)
			xw = zeros(3,1);
			xw(1) = this.x(1:3)'*xl + this.x(10);
			xw(2) = this.x(4:6)'*xl + this.x(11);
			xw(3) = this.x(7:9)'*xl + this.x(12);
		end

		%%
		function vw = transformVector(this,vl)
			vw = zeros(3,1);
			vw(1) = this.x(1:3)'*vl;
			vw(2) = this.x(4:6)'*vl;
			vw(3) = this.x(7:9)'*vl;
		end

		%%
		function xl = invTransformPoint(this,xw)
			A = reshape(this.x(1:9),3,3)';
			p = x(10:12);
			xl = A\(xw - p);
		end

		%%
		function vl = invTransformVector(this,vw)
			A = reshape(this.x(1:9),3,3)';
			vl = A\vw;
		end
		
		%%
		function computeInertiaConst(this)
			% Computes inertia for the shape
			I = this.shape.computeInertia(this.density);
            this.Wp = 1/I(4);
			R2A = 0.5*[
				-1  1  1
				 1 -1  1
				 1  1 -1
				];
			this.Wa = 1./(R2A*I(1:3));
		end

		%%
		function v = computePointVel(this,xl,k,ks,hs)
			xdot = this.computeVelocity(k,ks,hs);
			J = apbd.BodyAffine.jacobian(xl);
			v = J*xdot;
		end
		
        %%
		function stepBDF1(this,k,ks,hs,grav)
			xdot = this.computeVelocity(k,ks,hs);
            this.x0 = this.x;
			axdot = xdot(1:3);
			aydot = xdot(4:6);
			azdot = xdot(7:9);
			pdot = xdot(10:12);
			w = this.Wp;
			W = this.Wa;
			f = zeros(3,1);
			t = zeros(9,1); % affine torque
			tx = t(1:3);
			ty = t(4:6);
			tz = t(7:9);
			f = f + grav/w; % Gravity
			% Integrate velocities
			axdot = axdot + hs*(W.*tx);
			aydot = aydot + hs*(W.*ty);
			azdot = azdot + hs*(W.*tz);
			pdot  = pdot  + hs*(w *f);
			% Integrate positions
			this.x(1:3) = this.x(1:3) + hs*axdot;
			this.x(4:6) = this.x(4:6) + hs*aydot;
			this.x(7:9) = this.x(7:9) + hs*azdot;
			this.x(10:12) = this.x(10:12) + hs*pdot;
		end

		%%
		function [T,V] = computeEnergies(this,k,ks,hs,grav)
			xdot = this.computeVelocity(k,ks,hs);
			axdot = xdot(1:3);
			aydot = xdot(4:6);
			azdot = xdot(7:9);
			pdot = xdot(10:12);
			p = this.x(10:12);
			m = 1/this.Wp;
			I = 1./this.Wa;
			% Energy
			T = 0.5*(axdot'*(I.*axdot) + aydot'*(I.*aydot) + azdot'*(I.*azdot) + pdot'*m*pdot);
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
		end

		%%
		function [F,V] = draw(this)
			E = this.computeTransform();
			[F,V] = this.shape.draw(E,this.color,this.axisSize);
		end
	end

	methods (Static)
		%%
		function J = jacobian(xl)
			% Not used -- too expensive to form
			J = zeros(3,12);
			J(1,1:3) = xl';
			J(2,4:6) = xl';
			J(3,7:9) = xl';
			J(1:3,10:12) = eye(3);
		end

		%%
		function J = jacobianRot(xl)
			% Not used -- too expensive to form
			J = zeros(3,12);
			J(1,1:3) = xl';
			J(2,4:6) = xl';
			J(3,7:9) = xl';
		end
	end
end
