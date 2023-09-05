classdef ConstraintAffineJointRevolute < apbd.Constraint
	%ConstraintAffineJointRevolute Revolute joint between two affine bodies
	% This version uses two points, which are set 'height' apart

	properties
		affines
		xal1 % Point A for the joint in affine1's frame
		xal2 % Point A for the joint in affine2's frame
		xbl1 % Point B for the joint in affine1's frame
		xbl2 % Point B for the joint in affine2's frame

		% Drawing etc.
		radius
		height
	end

	methods
		%%
		function this = ConstraintAffineJointRevolute(affines,xls,als,height)
			this = this@apbd.Constraint(6,24);
			this.affines = affines;
			h = height/2;
			this.xal1 = xls{1} - h*als{1};
			this.xbl1 = xls{1} + h*als{1};
			if length(affines) == 2
				this.xal2 = xls{2} - h*als{2};
				this.xbl2 = xls{2} + h*als{2};
			end
			this.radius = 1;
			this.height = height;
		end

		%%
		function init(this)
			if length(this.affines) == 1
				% If only one body is provided, compute the world position
				% from the first local point.
				E = this.affines{1}.computeTransform();
				this.xal2 = E(1:3,:) * [this.xal1; 1];
				this.xbl2 = E(1:3,:) * [this.xbl1; 1];
			end
			% dC is constant, so we can precompute WdC and dCWdC
			for i = 1 : 6
				if length(this.affines) == 1
					affine1 = this.affines{1};
					W1 = affine1.W;
					xa1 = this.xal1;
					xb1 = this.xbl1;
					if i <= 3
						% Point A
						iR = 3*(i - 1) + (1:3); % indices for the rotational parts of q
						ip = 9 + i; % indices for the translationl part of q
						this.WdC(i,ip) = W1(4);
						this.WdC(i,iR) = [W1(1)*xa1(1), W1(2)*xa1(2), W1(3)*xa1(3)];
						this.dCWdC(i) = W1(4) + W1(1)*xa1(1)*xa1(1) + W1(2)*xa1(2)*xa1(2) + W1(3)*xa1(3)*xa1(3);
					else
						% Point B
						% We subtract 3 from i below because i=[4,5,6]
						% correspond to (x,y,z)
						iR = 3*((i - 3) - 1) + (1:3); % indices for the rotational parts of q
						ip = 9 + (i - 3); % indices for the translationl part of q
						this.WdC(i,ip) = W1(4);
						this.WdC(i,iR) = [W1(1)*xb1(1), W1(2)*xb1(2), W1(3)*xb1(3)];
						this.dCWdC(i) = W1(4) + W1(1)*xb1(1)*xb1(1) + W1(2)*xb1(2)*xb1(2) + W1(3)*xb1(3)*xb1(3);
					end
				else
					affine1 = this.affines{1};
					affine2 = this.affines{2};
					W1 = affine1.W;
					W2 = affine2.W;
					xa1 = this.xal1;
					xb1 = this.xbl1;
					xa2 = this.xal2;
					xb2 = this.xbl2;
					if i <= 3
						% Point A
						iR = 3*(i - 1) + (1:3);
						ip = 9 + i;
						this.WdC(i,ip) = W1(4);
						this.WdC(i,iR) = [W1(1)*xa1(1), W1(2)*xa1(2), W1(3)*xa1(3)];
						this.dCWdC(i) = W1(4) + W1(1)*xa1(1)*xa1(1) + W1(2)*xa1(2)*xa1(2) + W1(3)*xa1(3)*xa1(3);
						this.WdC(i,12+ip) = -W2(4);
						this.WdC(i,12+iR) = [-W2(1)*xa2(1), -W2(2)*xa2(2), -W2(3)*xa2(3)];
						this.dCWdC(i) = W2(4) + W2(1)*xa2(1)*xa2(1) + W2(2)*xa2(2)*xa2(2) + W2(3)*xa2(3)*xa2(3);
					else
						% Point B
						iR = 3*((i - 3) - 1) + (1:3);
						ip = 9 + (i - 3);
						this.WdC(i,ip) = W1(4);
						this.WdC(i,iR) = [W1(1)*xb1(1), W1(2)*xb1(2), W1(3)*xb1(3)];
						this.dCWdC(i) = W1(4) + W1(1)*xb1(1)*xb1(1) + W1(2)*xb1(2)*xb1(2) + W1(3)*xb1(3)*xb1(3);
						this.WdC(i,12+ip) = -W2(4);
						this.WdC(i,12+iR) = [-W2(1)*xb2(1), -W2(2)*xb2(2), -W2(3)*xb2(3)];
						this.dCWdC(i) = W2(4) + W2(1)*xb2(1)*xb2(1) + W2(2)*xb2(2)*xb2(2) + W2(3)*xb2(3)*xb2(3);
					end
				end
			end
		end

		%%
		function update(this,i)
			if length(this.affines) == 1
				affine1 = this.affines{1};
				if i <= 3
					% Point A
					iR = 3*(i - 1) + (1:3);
					ip = 9 + i;
					xaw1i = affine1.q(iR)'*this.xal1 + affine1.q(ip);
					this.C(i) = xaw1i - this.xal2(i);
				else
					% Point B
					iR = 3*((i - 3) - 1) + (1:3);
					ip = 9 + (i - 3);
					xbw1i = affine1.q(iR)'*this.xbl1 + affine1.q(ip);
					this.C(i) = xbw1i - this.xbl2(i - 3);
				end
			else
				affine1 = this.affines{1};
				affine2 = this.affines{2};
				if i <= 3
					% Point A
					iR = 3*(i - 1) + (1:3);
					ip = 9 + i;
					xaw1i = affine1.q(iR)'*this.xal1 + affine1.q(ip);
					xaw2i = affine2.q(iR)'*this.xal2 + affine2.q(ip);
					this.C(i) = xaw1i - xaw2i;
				else
					% Point B
					iR = 3*((i - 3) - 1) + (1:3);
					ip = 9 + (i - 3);
					xbw1i = affine1.q(iR)'*this.xbl1 + affine1.q(ip);
					xbw2i = affine2.q(iR)'*this.xbl2 + affine2.q(ip);
					this.C(i) = xbw1i - xbw2i;
				end
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
			T1(1:3,4) = 0.5*(this.xal1 + this.xbl1);
			R1 = eye(4);
			z = [0 0 1]';
			al1 = this.xbl1 - this.xal1;
			axis = cross(z,al1);
			axis = axis/norm(axis);
			angle = acos(z'*al1);
			R1(1:3,1:3) = se3.aaToMat(axis,angle);
			[X1,Y1,Z1] = se3.surfCylinder(E1*T1*R1,this.radius,this.height);
			surf(X1,Y1,Z1,'FaceColor',c1);
			if length(this.affines) == 2
				% Draw body 2
				affine2 = this.affines{2};
				E2 = affine2.computeTransform();
				c2 = affine2.color;
				T2 = eye(4);
				T2(1:3,4) = 0.5*(this.xal2 + this.xbl2);
				R2 = eye(4);
				al2 = this.xbl2 - this.xal2;
				axis = cross(z,al2);
				axis = axis/norm(axis);
				angle = acos(z'*al2);
				R2(1:3,1:3) = se3.aaToMat(axis,angle);
				[X2,Y2,Z2] = se3.surfCylinder(E2*T2*R2,0.5*this.radius,this.height);
				surf(X2,Y2,Z2,'FaceColor',c2);
			end
		end
	end
end
