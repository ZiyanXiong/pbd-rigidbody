classdef ConJointSph < apbd.ConJoint
	%ConJointSph Spherical joint between two rigid bodies
	%
	% Notation:
	% pl1 represents where the joint is wrt the body
	% (q1, p1) represent where the body is wrt the world
	% pw1 represents where the joint is wrt the world

	properties
		pl1 % Where the joint is in body1's frame
		pl2 % Where the joint is in body2's frame

		% Drawing etc.
		radius
	end

	methods
		%%
		function this = ConJointSph(bodies,xls)
			this = this@apbd.ConJoint(bodies,1);
			this.pl1 = xls{1};
			if length(bodies) == 2
				this.pl2 = xls{2};
			else
				this.pl2 = [];
			end
			this.radius = 1;
		end

		%%
		function init(this)
			if length(this.bodies) == 1
				% Constrain body1 to world, so set E2 to be where E1 is in
				% world frame
				this.pl2 = this.bodies{1}.transformPoint(this.pl1);
			end
		end

		%%
		function solve(this)
			if length(this.bodies) == 1
				rl1 = this.pl1; % in body1 coords
				rw2 = this.pl2; % already in world coords
				rw1 = this.bodies{1}.transformPoint(rl1); % in world coords
				dxw = rw1 - rw2;
				this.solvePos1(1,dxw,rl1);
			else
				rl1 = this.pl1; % in body1 coords
				rl2 = this.pl2; % in body2 coords
				rw1 = this.bodies{1}.transformPoint(rl1); % in world coords
				rw2 = this.bodies{2}.transformPoint(rl2); % in world coords
				dxw = rw1 - rw2;
				this.solvePos2(1,dxw,rl1,rl2);
			end
		end

		%%
		function draw(this)
			% Draw body 1
			body1 = this.bodies{1};
			E1 = body1.computeTransform();
			El1 = eye(4);
			El1(1:3,4) = this.pl1;
			c1 = body1.color;
			[X1,Y1,Z1] = se3.surfSphere(E1*El1,this.radius);
			surf(X1,Y1,Z1,'FaceColor',c1);
			if length(this.bodies) == 2
				% Draw body 2
				body2 = this.bodies{2};
				E2 = body2.computeTransform();
				El2 = eye(4);
				El2(1:3,4) = this.pl2;
				c2 = body2.color;
				[X2,Y2,Z2] = se3.surfSphere(E2*El2,0.5*this.radius);
				surf(X2,Y2,Z2,'FaceColor',c2);
			end
		end
	end
end
