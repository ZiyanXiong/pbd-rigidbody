classdef ConstraintAffineJointSpherical < apbd.Constraint
	%ConstraintAffineJointSpherical Spherical joint between two affine bodies

	properties
		affines
		xl1 % local position for the joint in body1's frame
		xl2 % local position for the joint in body2's frame

		% Drawing etc.
		radius
	end

	methods
		%%
		function this = ConstraintAffineJointSpherical(affines,xls)
			this = this@apbd.Constraint(3,24);
			this.affines = affines;
			this.xl1 = xls{1};
			if length(affines) == 2
				this.xl2 = xls{2};
			else
				this.xl2 = [];
			end
			this.radius = 1;
		end

		%%
		function init(this)
			% If only one body is provided, compute the world position from
			% the first local point.
			if length(this.affines) == 1
				E = this.affines{1}.computeTransform();
				this.xl2 = E(1:3,:) * [this.xl1; 1];
			end
			% dC is constant, so we can precompute WdC and dCWdC
			for i = 1 : 3
				ir = 3*(i - 1) + (1:3); % indices for the rotational parts of q
				ip = 9 + i; % indices for the translationl part of q
				if length(this.affines) == 1
					affine1 = this.affines{1};
					W1 = affine1.W;
					x1 = this.xl1;
					% We don't need the full Jacobian, so we form it
					% row by row
					%    dC = apbd.Body.jacobian(this.xl1);
					% is the same as
					%    dC(i,1) = 1;
					%    dC(i,4:6) = this.xl1';
					this.WdC(i,ip) = W1(4);
					this.WdC(i,ir) = [W1(1)*x1(1), W1(2)*x1(2), W1(3)*x1(3)];
					this.dCWdC(i) = W1(4) + W1(1)*x1(1)*x1(1) + W1(2)*x1(2)*x1(2) + W1(3)*x1(3)*x1(3);
				else
					affine1 = this.affines{1};
					affine2 = this.affines{2};
					W1 = affine1.W;
					W2 = affine2.W;
					x1 = this.xl1;
					x2 = this.xl2;
					this.WdC(i,ip) = W1(4);
					this.WdC(i,ir) = [W1(1)*x1(1), W1(2)*x1(2), W1(3)*x1(3)];
					dCWdCi = W1(4) + W1(1)*x1(1)*x1(1) + W1(2)*x1(2)*x1(2) + W1(3)*x1(3)*x1(3);
					this.WdC(i,12+ip) = -W2(4);
					this.WdC(i,12+ir) = [-W2(1)*x2(1), -W2(2)*x2(2), -W2(3)*x2(3)];
					dCWdCi = dCWdCi + W2(4) + W2(1)*x2(1)*x2(1) + W2(2)*x2(2)*x2(2) + W2(3)*x2(3)*x2(3);
					this.dCWdC(i) = dCWdCi;
				end
			end
		end

		%%
		function update(this,i)
			ir = 3*(i - 1) + (1:3); % indices for the rotational parts of q
			ip = 9 + i; % indices for the translationl part of q
			if length(this.affines) == 1
				affine1 = this.affines{1};
				% We don't need to form the whole 4x4 matrix, so we form it
				% row by row
				%    E = affine.computeTransform();
				%    xw = E1(1:3,:) * [xl; 1];
				% is the same as
				%    xw(1) = affine.q(1:3)'*xl + affine.q(10);
				xw1i = affine1.q(ir)'*this.xl1 + affine1.q(ip);
				this.C(i) = xw1i - this.xl2(i);
			else
				affine1 = this.affines{1};
				affine2 = this.affines{2};
				xw1i = affine1.q(ir)'*this.xl1 + affine1.q(ip);
				xw2i = affine2.q(ir)'*this.xl2 + affine2.q(ip);
				this.C(i) = xw1i - xw2i;
			end
		end

		%%
		function apply(this,i,dlambda)
			this.affines{1}.q = this.affines{1}.q + dlambda*this.WdC(i,1:12)';
			if length(this.affines) == 2
				this.affines{2}.q = this.affines{2}.q + dlambda*this.WdC(i,13:24)';
			end
		end

		%%
		function draw(this)
			% Draw body 1
			affine1 = this.affines{1};
			E1 = affine1.computeTransform();
			c1 = affine1.color;
			T1 = eye(4);
			T1(1:3,4) = this.xl1;
			[X1,Y1,Z1] = se3.surfSphere(E1*T1,this.radius);
			surf(X1,Y1,Z1,'FaceColor',c1);
			if length(this.affines) == 2
				% Draw body 2
				affine2 = this.affines{2};
				E2 = affine2.computeTransform();
				c2 = affine2.color;
				T2 = eye(4);
				T2(1:3,4) = this.xl2;
				[X2,Y2,Z2] = se3.surfSphere(E2*T2,0.5*this.radius);
				surf(X2,Y2,Z2,'FaceColor',c2);
			end
		end
	end
end
