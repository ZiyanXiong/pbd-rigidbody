classdef BodyRigid < apbd.Body
	%BodyRigid A rigid body with a quaternion for rotation
	
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
		function this = BodyRigid(shape,density)
			global CM; %#ok<GVMIS>

			this = this@apbd.Body(7);
			this.shape = shape;
			this.density = density;
			this.Mr = zeros(3,1);
			this.Mp = 0;
            this.w = zeros(3,1);
            this.v = zeros(3,1);
            this.deltaBody2Worldp = zeros(3,1);
            this.deltaBody2Worldq = zeros(4,1);
            this.deltaBody2Worldq(4) = 1;
            this.deltaLinDt = zeros(3,1);
            this.deltaAngDt = zeros(3,1);

			this.color = CM(mod(this.index-1,size(CM,1))+1,:);
			this.axisSize = 1;
		end
		
		%%
		function setInitTransform(this,E)
			this.xInit(1:4) = se3.matToQ(E(1:3,1:3));
			if this.xInit(4) < 0, this.xInit = -this.xInit; end
			this.xInit(5:7) = E(1:3,4);
		end

		%%
		function Einit = computeInitTransform(this)
			% This can be used only after calling setInitTransform()
			Einit = eye(4);
			Einit(1:3,1:3) = se3.qToMat(this.xInit(1:4));
            Einit(1:3,4) = this.xInit(5:7);
		end

		%%
		function setInitVelocity(this,phi)
			q = this.xInit(1:4);
			this.xdotInit(5:7) = se3.qRot(q,phi(4:6));
			this.xdotInit(1:4) = se3.wToQdot(q,phi(1:3));
            this.v = phi(4:6);
            this.w = phi(1:3);
			% Debug
			% Edot = E0*se3.brac(phi);
			% this.xdotInit(5:7) = Edot(1:3,4);
			% h = 1e-6;
			% R0 = E0(1:3,1:3);
			% R1 = R0*se3.exp(h*phi(1:3));
			% q0 = se3.matToQ(R0);
			% q1 = se3.matToQ(R1);
			% if q0(4) < 0, q0 = -q0; end
			% if q1(4) < 0, q1 = -q1; end
			% qdot_ = (q1 - q0)/h;
		end
		
        %%
        function updateStates(this, hs)
			%q = this.x0(1:4);
            %R = se3.qToMat(q);
            %invsqrtI = R * diag(sqrt(1./this.Mr)) * R';
            %angularMotionVel = invsqrtI * this.w;
            angularMotionVel = this.w;
            wNorm =  norm(angularMotionVel);
            if(wNorm>1e-9)
                halfWDt = 0.5 * wNorm * hs;
                dq = [angularMotionVel * sin(halfWDt)/ wNorm; 0];
                result = se3.qMul(dq, this.deltaBody2Worldq);
                result = result + this.deltaBody2Worldq * cos(halfWDt);
                this.deltaBody2Worldq = result / norm(result);
            end
            this.deltaBody2Worldp = this. deltaBody2Worldp + this.v * hs;

            this.deltaAngDt = this.deltaAngDt + this.w * hs;
            this.deltaLinDt = this.deltaLinDt + this.v * hs;

            this.x(1:4) = se3.qMul(this.deltaBody2Worldq, this.x0(1:4));
            this.x(5:7) = this.x0(5:7) + this.deltaBody2Worldp;
        end

        %%
        function integrateStates(this)
			%q = this.x0(1:4);
            %R = se3.qToMat(q);
            %invsqrtI = R * diag(sqrt(1./this.Mr)) * R';
            %this.w = invsqrtI * this.w;
            this.x(1:4) = se3.qMul(this.deltaBody2Worldq, this.x0(1:4));
            this.x(5:7) = this.x0(5:7) + this.deltaBody2Worldp;
            this.deltaBody2Worldp = zeros(3,1);
            this.deltaBody2Worldq = zeros(4,1);
            this.deltaBody2Worldq(4) = 1;
            
            this.deltaLinDt = zeros(3,1);
            this.deltaAngDt = zeros(3,1);

            %Clear contact information
            this.layer = 99;
            this.neighbors = [];
        end

		%%
		function E = computeTransform(this)
			E = eye(4);
			E(1:3,1:3) = se3.qToMat(this.x(1:4));
            E(1:3,4) = this.x(5:7);
		end

		%%
		function xw = transformPoint(this,xl)
			q = this.x(1:4);
			p = this.x(5:7);
			xw = se3.qRot(q,xl) + p;
        end

		%%
		function xw = transformPointPrev(this,xl)
			q = this.x0(1:4);
			p = this.x0(5:7);
			xw = se3.qRot(q,xl) + p;
		end

		%%
		function vw = transformVector(this,vl)
			q = this.x(1:4);
			vw = se3.qRot(q,vl);
		end

		%%
		function xl = invTransformPoint(this,xw)
			q = this.x(1:4);
			p = this.x(5:7);
			xl = se3.qRotInv(q,xw - p);
		end

		%%
		function vl = invTransformVector(this,vw)
			q = this.x(1:4);
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
		function v = computePointVel(this,xl)
            %{
			xdot = this.computeVelocity(k,ks,hs);
			qdot = xdot(1:4);
			pdot = xdot(5:7); % in world coords
			q = this.x(1:4);
			w = se3.qdotToW(q,qdot); % angular velocity in body coords
			% v = R*cross(w,xl) + pdot
			v = se3.qRot(q,se3.cross(w,xl)) + pdot;
            %}
            q = this.x(1:4);
            rw = se3.qRot(q,xl);
            v = se3.cross(this.w,rw) + this.v;
		end
		
        %%
		function stepBDF1(this,h,grav)
            this.x0 = this.x;
			v = this.v; % pdot
			q = this.x(1:4);
			p = this.x(5:7);
            R = se3.qToMat(q);
            %R = eye(3);
			w = this.w; % angular velocity in body coords
			f = zeros(3,1); % translational force in world space
			t = zeros(3,1); % angular torque in body space
			m = this.Mp; % scalar mass
			I = R * diag(this.Mr) * R'; % inertia in world space;
			Iw = I*w; % angular momentum in body space
			f = f + m*grav; % Gravity
			t = t + se3.cross(Iw,w); % Coriolis
			% Integrate velocities
			w = w + h*(I\t);
			v = v + h*(m \f);
			%qdot = se3.wToQdot(q,w);
			% Integrate positions
			%q = q + hs*qdot;
			%p = p + hs*v;
			%q = q/norm(q);
			%this.x(1:4) = q;
			%this.x(5:7) = p;
            %sqrtInertia = R * diag(sqrt(this.Mr)) * R';
            %this.w = sqrtInertia * w;
            this.w = w;
            this.v = v;
		end

		%%
		function [T,V] = computeEnergies(this,k,ks,hs,grav)
			%xdot = this.computeVelocity(k,ks,hs);
			%qdot = xdot(1:4);
			%v = xdot(5:7); % pdot
            v = this.v; % pdot
            w = this.w; % angular velocity in body coords
			q = this.x(1:4);
			p = this.x(5:7);
			%w = se3.qdotToW(q,qdot); % angular velocity in body coords
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
