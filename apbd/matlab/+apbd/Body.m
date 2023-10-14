classdef (Abstract) Body < handle
	%Body An abstract body
	
	%%
	properties
		n        % DOF count
		xInit    % Initial position
		xdotInit % Initial velocity
        x        % position
        x0       % Previous position
        x1       % Position after BDF1 Update
		dxJacobi % Jacobi updates
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
			this.dxJacobi = zeros(n,1);
			this.collide = false;
			this.mu = 0;

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
			this.x = this.x + this.dxJacobi;
			this.dxJacobi = zeros(this.n,1);
            %{
			E = this.computeTransform();
            if ~(det(E(1:3,1:3))==1)
                print("Not Rigid")
                q = apbd.BodyRigid2d.unproj(this.x);
                qNormalized = q ./ norm(q,2);
                this.x(1:2) = qNormalized(3:4);
                E = this.computeTransform();
            end
            %}
            
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
