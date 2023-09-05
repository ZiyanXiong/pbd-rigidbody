classdef ConstraintJointRevolute5d < apbd.Constraint
	%ConstraintJointRevolute Constraint for revolute joint
	% This version uses a point and a direction.

	properties
		xl1 % local position for the joint in body1's frame
		xl2 % local position for the joint in body2's frame
		al1 % local z axis for the joint in body1's frame
		al2 % local xy axes for the joint in body2's frame

		% Drawing etc.
		radius
		height
	end

	methods
		%%
		function this = ConstraintJointRevolute5d(bodies,xls,als)
			this = this@apbd.Constraint(5,bodies);
			this.xl1 = xls{1};
			this.al1 = als{1}/norm(als{1});
			if length(bodies) == 2
				this.xl2 = xls{2};
				% We will compute xy from this z in init()
				this.al2 = als{2}/norm(als{2});
			else
				this.xl2 = [];
				this.al2 = [];
			end
			this.radius = 1;
			this.height = 1;
		end

		%%
		function init(this)
			% If only one body is provided, compute the world position from
			% the first local point.
			if length(this.bodies) == 1
				E = this.bodies{1}.computeTransform();
				this.xl2 = E(1:3,:) * [this.xl1; 1];
				this.al2 = E(1:3,1:3) * this.al1;
			end
			% Compute [x2 y2] for body2 so that [x2 y2 z2] form an
			% orthogonal basis.
			z2 = this.al2;
			x2 = se3.cross([0 1 0]',z2);
			if norm(x2) < 1e-6
				% x2 is zero. Generate new x2.
				if abs(z2) < 1e-6
					tmp = [0 0 1]';
				else
					tmp = [1 0 0]';
				end
				x2 = se3.cross(tmp,z2);
			end
			x2 = x2/norm(x2);
			y2 = se3.cross(z2,x2);
			this.al2 = [x2 y2];
		end

		%%
		function update(this,i)
			ii = 3+3*(i-1)+(1:3); % indices for the rotational part of q
			if length(this.bodies) == 1
				body1 = this.bodies{1};
				if i <= 3
					% Position
					xw1i = body1.q(ii)'*this.xl1 + body1.q(i);
					this.C(i) = xw1i - this.xl2(i);
					this.dC(i,i) = 1;
					this.dC(i,ii) = this.xl1';
				else
					% Orientation
					R1 = reshape(body1.q(4:12),3,3)';
					az1 = R1 * this.al1;
					Jz1 = apbd.Body.jacobianRot(az1);
					if i == 4
						ax2 = this.al2(1:3,1);
						this.C(4)    = ax2' * az1;
						this.dC(4,:) = ax2' * Jz1;
					else
						ay2 = this.al2(1:3,2);
						this.C(5)    = ay2' * az1;
						this.dC(5,:) = ay2' * Jz1;
					end
				end
			else
				body1 = this.bodies{1};
				body2 = this.bodies{2};
				if i <= 3
					% Position
					body1 = this.bodies{1};
					body2 = this.bodies{2};
					xw1i = body1.q(ii)'*this.xl1 + body1.q(i);
					xw2i = body2.q(ii)'*this.xl2 + body2.q(i);
					this.C(i) = xw1i - xw2i;
					this.dC(i,i) = 1;
					this.dC(i,ii) = this.xl1';
					this.dC(i,12+i) = -1;
					this.dC(i,12+ii) = -this.xl2';
				else
					% Orientation
					R1 = reshape(body1.q(4:12),3,3)';
					R2 = reshape(body2.q(4:12),3,3)';
					az1 = R1 * this.al1;
					Jz1 = apbd.Body.jacobianRot(az1);
					if i == 4
						ax2 = R2 * this.al2(1:3,1);
						this.C(4) = ax2' * az1;
						Jx2 = apbd.Body.jacobianRot(this.al2(1:3,1));
						this.dC(4,:) = [ax2' * Jz1, az1' * Jx2];
					else
						ay2 = R2 * this.al2(1:3,2);
						this.C(5) = ay2' * az1;
						Jy2 = apbd.Body.jacobianRot(this.al2(1:3,2));
						this.dC(5,:) = [ay2' * Jz1, az1' * Jy2];
					end
				end
			end
		end

		%%
		function draw(this)
			% Draw body 1
			body1 = this.bodies{1};
			E1 = body1.computeTransform();
			c1 = body1.color;
			T1 = eye(4);
			T1(1:3,4) = this.xl1;
			R1 = eye(4);
			z = [0 0 1]';
			axis = cross(z,this.al1);
			axis = axis/norm(axis);
			angle = acos(z'*this.al1);
			R1(1:3,1:3) = se3.aaToMat(axis,angle);
			[X1,Y1,Z1] = se3.surfCylinder(E1*T1*R1,this.radius,this.height);
			surf(X1,Y1,Z1,'FaceColor',c1);
			if length(this.bodies) == 2
				% Draw body 2
				body2 = this.bodies{2};
				E2 = body2.computeTransform();
				c2 = body2.color;
				T2 = eye(4);
				T2(1:3,4) = this.xl2;
				R2 = eye(4);
				z2 = se3.cross(this.al2(1:3,1),this.al2(1:3,2));
				axis = cross(z,z2);
				axis = axis/norm(axis);
				angle = acos(z'*z2);
				R2(1:3,1:3) = se3.aaToMat(axis,angle);
				[X2,Y2,Z2] = se3.surfCylinder(E2*T2*R2,0.5*this.radius,this.height);
				surf(X2,Y2,Z2,'FaceColor',c2);
			end
		end
	end
end
