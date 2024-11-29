classdef (Abstract) Body < handle
	%Body An abstract body
	
	%%
	properties
		n        % DOF count
		xInit    % Initial position
		xdotInit % Initial velocity
        x        % position
        x0       % Previous position
		dxJacobi % Jacobi updates
        dphiJacobi % Jacobi updates
		collide  % Collision on/off
		mu       % Coefficient of friction

        v        % Linear Velocity
        w        % Angular Velocity
        deltaLinDt % Change of linear motion
        deltaAngDt % Change of angular motion
        deltaBody2Worldp % Change of linear motion
        deltaBody2Worldq % Change of angular motion

		% Drawing etc.
		index    % Body index
        layer    % Which contact layer this body is
        neighbors % Indices for neighobor bodies
    end

	%%
	methods
		%%
		function this = Body(n)
			this.n = n;
			this.xInit = zeros(n,1);
			this.xdotInit = zeros(n,1);
			this.x = zeros(n,1);
			this.x0 = zeros(n,1);
			this.dxJacobi = zeros(n,1);
            this.dphiJacobi = zeros(n-1,1);
			this.collide = false;
			this.mu = 0;

			this.index = apbd.Model.countB();
            this.layer = 99;
            this.neighbors = [];
		end

		%%
		function init(this)
			this.computeInertiaConst();
			this.x = this.xInit;
			this.x0 = this.x;
		end

		%%
		function clearJacobi(this)
			this.dxJacobi = zeros(this.n,1);
		end

		%%
		function applyJacobi(this)
			this.x = this.x + this.dxJacobi;
            if length(this.x) == 7
                q = this.x(1:4);
                qNormalized = q ./ norm(q,2);
                this.x(1:4) = qNormalized;
            end
			this.dxJacobi = zeros(this.n,1);
        end

		%%
		function applyVelJacobi(this)
            if length(this.x) == 7
                this.v = this.v + this.dphiJacobi(1:3);
                this.w = this.w + this.dphiJacobi(4:6);
            end
			this.dphiJacobi = zeros(this.n-1,1);
		end

		%%
		function xdot = computeVelocity(this,k,ks,hs)
            %{
			% Compute velocity from position
			if k == 0 && ks == 0
				xdot = this.xdotInit;
			else
				xdot = (this.x - this.x0)/hs;
            end
			qdot = xdot(1:4);
			this.v = xdot(5:7); % in world coords
			q = this.x(1:4);
			this.w = se3.qdotToW(q,qdot); % angular velocity in body coords
            %}
            
            %{
			% Multiply method
			if k == 0 && ks == 0
				xdot = this.xdotInit;
			else
				xdot = (this.x - this.x0)/hs;
                xdot(1:4) = se3.qMulInv(this.x(1:4), this.x0(1:4));
            end
            %}

			xdot = (this.x - this.x0)/hs;
            xdot(1:4) = se3.qMulInv(this.x(1:4), this.x0(1:4));
			dq = xdot(1:4);
			this.v = xdot(5:7); % in world coords
			this.w = 2* dq(1:3) / hs;
            if dq(4) < 0
                this.w = -this.w;
            end
            
		end
            

	end

	methods (Abstract)
		%% Compute constant inertia if applicable
		computeInertiaConst(this);

		%%
		stepBDF1(this,k,ks,hs,grav);

		%%
		[T,V] = computeEnergies(this,k,ks,hs,grav);

		%%
		flag = broadphaseGround(this,Eg);

		%%
		cdata = narrowphaseGround(this,Eg);

		%%
		flag = broadphaseRigid(this,that);

		%%
		cdata = narrowphaseRigid(this,that);

		%%
		[F,V] = draw(this);
	end
end
