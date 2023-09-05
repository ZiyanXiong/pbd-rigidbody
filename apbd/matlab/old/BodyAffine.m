classdef BodyAffine < apbd.Body
	%Affine An affine body
	
	%%
	properties
		shape    % Associated shape
		density  % Mass/volume
		M        % Affine inertia (top 3 for rotation, bottom 1 for translation)
        W        % Inverse affine inertia (top 3 for rotation, bottom 1 for translation)

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
			this.M = zeros(4,1); % top 3 for rotation, bottom 1 for translation
			this.W = zeros(4,1); % top 3 for rotation, bottom 1 for translation

			this.color = CM(mod(this.index-1,size(CM,1))+1,:);
			this.axisSize = 1;
		end
		
		%%
		function setInitTransform(this,E)
			% Sets the transform of this Affine
            % Note: whenever we use reshape(), don't forget that we want q to be
            % row-wise, so we need the transpose.
            this.qInit = [reshape(E(1:3,1:3)',9,1); E(1:3,4)];
		end

		%%
		function setInitVelocity(this,phi)
			E = eye(4);
			E(1:3,1:3) = reshape(this.qInit(1:9),3,3)';
			E(1:3,4) = this.qInit(10:12);
			Edot = E*se3.brac(phi);
			this.qdotInit(1:9) = reshape(Edot(1:3,1:3)',9,1);
			this.qdotInit(10:12) = Edot(1:3,4);
		end
		
		%%
		function E = computeTransform(this)
			E = eye(4);
            E(1:3,1:3) = reshape(this.q(1:9),3,3)'; % row-by-row
            E(1:3,4) = this.q(10:12);
		end
		
		%%
		function computeInertiaConst(this)
			% Computes inertia for the shape
			I = this.shape.computeInertia(this.density);
            % Rotational inertia
            % This matrix transforms from rigid to affine.
            R2A = [
	            -1  1  1
	             1 -1  1
	             1  1 -1
	            ]/2;
			this.M(1:3,1) = R2A*I(1:3);
            % Translational inertia
            this.M(4,1) = I(4);
			% Invert
			this.W = 1./this.M;
		end

		%%
		function f = computeInertiaQVV(this,qdot) %#ok<INUSD>
			% Inertia is constant, so the quadratic velocity vector is zero
			f = zeros(this.n,1);
		end

        %%
		function a = computeUnconAcc(this,grav,f)
			f(10:12) = f(10:12) + this.M(4)*grav;
			a = zeros(12,1);
			a( 1: 3) = this.W(1:3).*f( 1: 3);
			a( 4: 6) = this.W(1:3).*f( 4: 6);
			a( 7: 9) = this.W(1:3).*f( 7: 9);
            a(10:12) = this.W(4)*f(10:12);
		end

		%%
		function postStep(this) %#ok<MANU>
			% Do nothing
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
