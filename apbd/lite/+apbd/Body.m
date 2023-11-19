classdef (Abstract) Body < handle
	%Body An abstract body
	
	%%
	properties
		n        % DOF count
		xInit    % Initial position
		xdotInit % Initial velocity
        x        % Projected position
        x0       % Previous position
        x1       % Positions updated in solvers
        x1_0     % Position after BDF1 Update
		dxJacobi % Jacobi updates
        dxJacobiShock % Jacobi updates after shock propagation
        shockParentConIndex % Constraints for shock propagation
        conIndex % Constraints related to this body
        layer    % How far from the ground
		collide  % Collision on/off
		mu       % Coefficient of friction

		% Drawing etc.
		index    % Body index
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
            this.x1 = zeros(n,1);
            this.x1_0 = zeros(n,1);
			this.dxJacobi = zeros(n,1);
            this.dxJacobiShock = zeros(n,1);
			this.collide = false;
			this.mu = 0;
            this.shockParentConIndex = {};
            this.conIndex = [];
            this.layer = 99;

			this.index = apbd.Model.countB();
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
			this.x1 = this.x1 + this.dxJacobi;
			this.dxJacobi = zeros(this.n,1);
            this.regularize();
        end

		%%
		function applyJacobiShock(this)
			this.x1 = this.x1 + this.dxJacobiShock;
			this.dxJacobiShock = zeros(this.n,1);
            this.regularize();
        end

        %%
        function regularize(this)
            this.x = this.x1;
            if length(this.x) == 4
                q = apbd.BodyRigid2d.unproj(this.x);
                qNormalized = q ./ norm(q,2);
                this.x(1:2) = qNormalized(3:4);
            elseif length(this.x) == 7
                q = this.x(1:4);
                qNormalized = q ./ norm(q,2);
                this.x(1:4) = qNormalized;
            elseif length(this.x) == 12
                A = reshape(this.x(1:9),3,3)';
                [U,~,V] = svd(A);
                A = U * V';
                this.x(1:9) = reshape(A',9,1);
            end
        end

		%%
		function xdot = computeVelocity(this,k,ks,hs)
			% Compute velocity from position
			if k == 0 && ks == 0
				xdot = this.xdotInit;
			else
				xdot = (this.x - this.x0)/hs;
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
