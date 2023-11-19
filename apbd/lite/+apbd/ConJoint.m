classdef (Abstract) ConJoint < apbd.ConBase
	%ConJoint Joint between two rigid bodies
	%
	% Notation:
	% pl1 represents where the joint is wrt the body
	% (q1, p1) represent where the body is wrt the world
	% pw1 represents where the joint is wrt the world

	properties
		bodies
        shockProp % Enable shock propagation or not
	end

	methods
		%%
		function this = ConJoint(bodies,m)
			this = this@apbd.ConBase(m);
			this.bodies = bodies;
            this.shockProp = false;
		end

		%%
		function [dlambda,dq1] = solveRot1(this,i,dqw,apply)
			% [Muller et al. 2020] Section 3.3.2 for one body
			if nargin < 4
				apply = true;
			end
			dlambda = 0;
			dq1 = zeros(4,1);
			this.C(i) = norm(dqw);
			if this.C(i) > 1e-6
				body1 = this.bodies{1};
				m1 = body1.Mp;
				I1 = body1.Mr;
				q1 = body1.x(1:4);
				nw = dqw/this.C(1);
				nl1 = se3.qRotInv(q1,nw);
				w1 = (1/m1) + nl1'*(I1.\nl1);
				numerator = -this.C(i);
				denominator = w1;
				dlambda = numerator/denominator;
				this.lambda(i) = this.lambda(i) + dlambda;
				% Quaternion update
                q1 = body1.x1_0(1:4);
				dpw = dlambda*nw;
				dpl1 = se3.qRotInv(q1,dpw);
				qtmp1 = [se3.qRot(q1,I1.\dpl1); 0];
				dq1 = 0.5*se3.qMul(qtmp1,q1);
				if apply
					body1.x1(1:4) = body1.x1(1:4) + dq1;
                    body1.regularize();
				end
			end
		end

		%%
		function [dlambda,dq1,dq2] = solveRot2(this,i,dqw,apply)
			% [Muller et al. 2020] Section 3.3.2 for two bodies
			if nargin < 4
				apply = true;
			end
			dlambda = 0;
			dq1 = zeros(4,1);
			dq2 = zeros(4,1);
			this.C(i) = norm(dqw);
			if this.C(i) > 1e-6
				body1 = this.bodies{1};
				body2 = this.bodies{2};
				m1 = body1.Mp;
				m2 = body2.Mp;
				I1 = body1.Mr;
				I2 = body2.Mr;
				q1 = body1.x(1:4);
				q2 = body2.x(1:4);
				nw = dqw/this.C(i);
				nl1 = se3.qRotInv(q1,nw);
				nl2 = se3.qRotInv(q2,nw);
				w1 = (1/m1) + nl1'*(I1.\nl1);
				w2 = (1/m2) + nl2'*(I2.\nl2);
				numerator = -this.C(i);
				denominator = w1 + w2;
				dlambda = numerator/denominator;
				this.lambda(i) = this.lambda(i) + dlambda;
				% Quaternion update
                q1 = body1.x1_0(1:4);
                q2 = body2.x1_0(1:4);
				dpw = dlambda*nw;
				dpl1 = se3.qRotInv(q1,dpw);
				dpl2 = se3.qRotInv(q2,dpw);
				qtmp1 = [se3.qRot(q1,I1.\dpl1); 0];
				qtmp2 = [se3.qRot(q2,I2.\dpl2); 0];
				dq1 =  0.5*se3.qMul(qtmp1,q1);
				dq2 = -0.5*se3.qMul(qtmp2,q2);
				if apply
                    if this.shockProp
		                body1.dxJacobiShock(1:4) = body1.dxJacobiShock(1:4) + dq1;
                    else
		                body1.x1(1:4) = body1.x1(1:4) + dq1;
                    end
					body1.x1(1:4) = body1.x1(1:4) + dq1;
					body2.x1(1:4) = body2.x1(1:4) + dq2;
                    body1.regularize();
                    body2.regularize();
				end
			end
		end

		%%
		function [dlambda,dq1,dp1] = solvePos1(this,i,dxw,rl1,apply)
			% [Muller et al. 2020] Section 3.3.1 for one body
			if nargin < 5
				apply = true;
			end
			dlambda = 0;
			dq1 = zeros(4,1);
			dp1 = zeros(3,1);
			this.C(i) = norm(dxw);
			if this.C(i) > 1e-6
				body1 = this.bodies{1};
				m1 = body1.Mp;
				I1 = body1.Mr;
				q1 = body1.x(1:4);
				nw = dxw/this.C(i);
				nl1 = se3.qRotInv(q1,nw);
				rnl1 = se3.cross(rl1,nl1);
				w1 = (1/m1) + rnl1'*(I1.\rnl1);
				numerator = -this.C(i);
				denominator = w1;
				dlambda = numerator/denominator;
				this.lambda(i) = this.lambda(i) + dlambda;
				% Position update
				dpw = dlambda*nw;
				dp1 = dpw/m1;
				% Quaternion update
                q1 = body1.x1_0(1:4);
				dpl1 = se3.qRotInv(q1,dpw);
				qtmp1 = [se3.qRot(q1,I1.\se3.cross(rl1,dpl1)); 0];
				dq1 = 0.5*se3.qMul(qtmp1,q1);
				if apply
					body1.x1(1:4) = body1.x1(1:4) + dq1;
                    body1.x1(5:7) = body1.x1(5:7) + dp1;
                    body1.regularize();
				end
			end
		end

		%%
		function [dlambda,dq1,dp1,dq2,dp2] = solvePos2(this,i,dxw,rl1,rl2,apply)
			% [Muller et al. 2020] Section 3.3.1 for two bodies
			if nargin < 6
				apply = true;
			end
			dlambda = 0;
			dq1 = zeros(4,1);
			dp1 = zeros(3,1);
			dq2 = zeros(4,1);
			dp2 = zeros(3,1);
			this.C(i) = norm(dxw);
			if this.C(i) > 1e-6
				body1 = this.bodies{1};
				body2 = this.bodies{2};
				m1 = body1.Mp;
				m2 = body2.Mp;
				I1 = body1.Mr;
				I2 = body2.Mr;
				q1 = body1.x(1:4);
				q2 = body2.x(1:4);
				nw = dxw/this.C(i);
				nl1 = se3.qRotInv(q1,nw);
				nl2 = se3.qRotInv(q2,nw);
				rnl1 = se3.cross(rl1,nl1);
				rnl2 = se3.cross(rl2,nl2);
				w1 = (1/m1) + rnl1'*(I1.\rnl1);
				w2 = (1/m2) + rnl2'*(I2.\rnl2);
				numerator = -this.C(i);
				denominator = w1 + w2;
				dlambda = numerator/denominator;
				this.lambda(i) = this.lambda(i) + dlambda;
				% Position update
				dpw = dlambda*nw;
				dp1 =  dpw/m1;
				dp2 = -dpw/m2;
				% Quaternion update
                q1 = body1.x1_0(1:4);
                q2 = body2.x1_0(1:4);
				dpl1 = se3.qRotInv(q1,dpw);
				dpl2 = se3.qRotInv(q2,dpw);
				qtmp1 = [se3.qRot(q1,I1.\se3.cross(rl1,dpl1)); 0];
				qtmp2 = [se3.qRot(q2,I2.\se3.cross(rl2,dpl2)); 0];
				dq1 =  0.5*se3.qMul(qtmp1,q1);
				dq2 = -0.5*se3.qMul(qtmp2,q2);
				if apply
                    if this.shockProp
		                body1.dxJacobiShock(1:4) = body1.dxJacobiShock(1:4) + dq1;
                        body1.dxJacobiShock(5:7) = body1.dxJacobiShock(5:7) + dp1;
                    else
					    body1.x1(1:4) = body1.x1(1:4) + dq1;
					    body1.x1(5:7) = body1.x1(5:7) + dp1;
                    end
					body2.x1(1:4) = body2.x1(1:4) + dq2;
				    body2.x1(5:7) = body2.x1(5:7) + dp2;
                    body1.regularize();
                    body2.regularize();
				end
			end
		end
	end
end
