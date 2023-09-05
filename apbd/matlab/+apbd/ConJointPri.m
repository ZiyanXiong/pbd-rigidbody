classdef ConJointPri < apbd.ConJoint
	%ConJointPri Prismatic joint between two rigid bodies
	%
	% Notation:
	% (ql1, pl1) represent where the joint is wrt the body
	% (q1,  p1)  represent where the body is wrt the world
	% (qw1, pw1) represent where the joint is wrt the world

	properties
		ql1 % Joint orientation wrt body1's frame
		pl1 % Joint position wrt body1's frame
		ql2 % Joint orientation wrt body2's frame
		pl2 % Joint position wrt body2's frame

		% Drawing etc.
		width
		height
	end

	methods
		%%
		function this = ConJointPri(bodies,Es)
			this = this@apbd.ConJoint(bodies,2);
			E1 = Es{1};
			this.ql1 = se3.matToQ(E1(1:3,1:3));
			if this.ql1(4) < 0, this.ql1 = -this.ql1; end
			this.pl1 = E1(1:3,4);
			if length(bodies) == 2
				E2 = Es{2};
				this.ql2 = se3.matToQ(E2(1:3,1:3));
				if this.ql2(4) < 0, this.ql2 = -this.ql2; end
				this.pl2 = E2(1:3,4);
			else
				this.ql2 = [];
				this.pl2 = [];
			end
		end

		%%
		function init(this)
			if length(this.bodies) == 1
				% Constrain body1 to world, so set E2 to be where E1 is in
				% world frame
				E1 = this.bodies{1}.computeTransform();
				El1 = eye(4);
				El1(1:3,1:3) = se3.qToMat(this.ql1);
				El1(1:3,4) = this.pl1;
				E2 = E1*El1;
				this.ql2 = se3.matToQ(E2(1:3,1:3));
				if this.ql2(4) < 0, this.ql2 = -this.ql2; end
				this.pl2 = E2(1:3,4);
			end
		end

		%%
		function solve(this)
			if length(this.bodies) == 1
				body1 = this.bodies{1};
				q1 = body1.x(1:4);
				p1 = body1.x(5:7);
				% Rotation
				qw1 = se3.qMul(q1,this.ql1);
				dqw = se3.qMulInv(qw1,this.ql2);
				dqw = dqw(1:3);
				this.solveRot1(1,dqw);
				% Position
				rl1 = this.pl1;
				rw1 = se3.qRot(q1,rl1) + p1;
				drw = rw1 - this.pl2; % pl2 is in world coords
				aw2 = se3.qToMatZ(this.ql2); % unconstrained direction in world coords
				a = aw2'*drw;
				drw = drw - a*aw2;
				this.solvePos1(2,drw,rl1);
			else
				body1 = this.bodies{1};
				body2 = this.bodies{2};
				q1 = body1.x(1:4);
				q2 = body2.x(1:4);
				p1 = body1.x(5:7);
				p2 = body2.x(5:7);
				% Rotation
				qw1 = se3.qMul(q1,this.ql1);
				qw2 = se3.qMul(q2,this.ql2);
				dqw = se3.qMulInv(qw1,qw2);
				dqw = dqw(1:3);
				this.solveRot2(1,dqw);
				% Position
				rl1 = this.pl1;
				rl2 = this.pl2;
				rw1 = se3.qRot(q1,rl1) + p1;
				rw2 = se3.qRot(q2,rl2) + p2;
				drw = rw1 - rw2;
				qw1 = se3.qMul(q1,this.ql1);
				aw1 = se3.qToMatZ(qw1);
				a = aw1'*drw;
				drw = drw - a*aw1;
				this.solvePos2(2,drw,rl1,rl2);
			end
		end

		%%
		function draw(this)
			% Draw body 1
			body1 = this.bodies{1};
			E1 = body1.computeTransform();
			El1 = eye(4);
			El1(1:3,1:3) = se3.qToMat(this.ql1);
			El1(1:3,4) = this.pl1;
			Ew1 = E1*El1;
			[F,V] = se3.patchCuboid(Ew1,[this.width this.width this.height]);
			patch('Faces',F,'Vertices',V,'FaceColor',body1.color);
			if length(this.bodies) == 2
				% Draw body 2
				body2 = this.bodies{2};
				E2 = body2.computeTransform();
				El2 = eye(4);
				El2(1:3,1:3) = se3.qToMat(this.ql2);
				El2(1:3,4) = this.pl2;
				Ew2 = E2*El2;
				[F,V] = se3.patchCuboid(Ew2,[this.width this.width this.height]);
				patch('Faces',F,'Vertices',V,'FaceColor',body2.color);
			end
		end
	end
end
