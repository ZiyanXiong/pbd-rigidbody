classdef ConCollRigidRigid < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
		body1
		body2
		x1 % Position wrt body 1 (3x1)
		x2 % Position wrt body 2 (3x1)
	end

	methods
		%%
		function this = ConCollRigidRigid(body1,body2)
			this = this@apbd.ConColl();
			this.body1 = body1;
			this.body2 = body2;
            this.bodies = {body1;body2};
			this.x1 = zeros(3,1);
			this.x2 = zeros(3,1);
		end

		%%
		function init(this) %#ok<MANU>
			% Do nothing
		end

		%%
		function setData(this,cdata)
			this.d = cdata.d;
			this.nw = cdata.nw;
			this.x1 = cdata.x1;
			this.x2 = cdata.x2;
        end

        %%
        function applyJacobi(this)
            this.body1.applyJacobi();
            this.body2.applyJacobi();
        end

        %%
        function swapBody(this)
            this.nw = -this.nw;
            temp = this.body2;
            this.body2 = this.body1;
            this.body1 = temp;

            this.bodies{1} = this.body1;
            this.bodies{2} = this.body2;

            temp = this.x2;
            this.x2 = this.x1;
            this.x1 = temp;
        end

        %%
        function [dq1,dp1,dq2,dp2] = computeDx(this, dlambda, nw)
			m1 = this.body1.Mp;
			I1 = this.body1.Mr;
			m2 = this.body2.Mp;
			I2 = this.body2.Mr;
			% Position update
			dpw = dlambda*nw;
			dp1 =  dpw/m1;
			dp2 = -dpw/m2;
			% Quaternion update
            q1 = this.body1.x1_0(1:4);
            q2 = this.body2.x1_0(1:4);
			dpl1 = se3.qRotInv(q1,dpw);
			dpl2 = se3.qRotInv(q2,dpw);
			qtmp1 = [se3.qRot(q1,I1.\se3.cross(this.x1,dpl1)); 0];
			qtmp2 = [se3.qRot(q2,I2.\se3.cross(this.x2,dpl2)); 0];
            %dq1 = se3.qMul(sin(0.5*qtmp1),q1);
            %dq2 = se3.qMul(sin(-0.5*qtmp2),q2);
		    dq1 =  0.5*se3.qMul(qtmp1,q1);
			dq2 = -0.5*se3.qMul(qtmp2,q2);
        end

		%%
		function solveNorPos(this, hs)
			v1w = this.body1.computePointVel(this.x1, hs);
			v2w = this.body2.computePointVel(this.x2, hs);

            % 0.1 controls how much penetration solved at each step, 1e-3
            % is the penetration we hope to keep. These two parameters can
            % be adjusted according to scene configuration
			v = hs * (v1w - v2w) - 0.1 * (this.d + 1e-3) * this.nw;
            vNorm = norm(v);
            vNormalized = v ./ vNorm;
            [tx,ty] = apbd.ConColl.generateTangents(this.nw);
            vNormalizedContactFrame = [-this.nw'; tx' ; ty'] * vNormalized;

            dlambda = this.solvePosDir2(norm(v), vNormalized);
            this.C = vNorm * vNormalizedContactFrame;

            dlambdaNor = dlambda * vNormalizedContactFrame(1);
            lambdaNor = this.lambda(1) + dlambdaNor;
            if lambdaNor < 0
                dlambdaNor  = - this.lambda(1);
            end
            this.lambda(1) = this.lambda(1) + dlambdaNor;
			mu1 = this.body1.mu;
			mu2 = this.body2.mu;
			mu = 0.5*(mu1 + mu2);
            dlambdaTan = [0;0];
            if mu > 0
                dlambdaTx = dlambda * vNormalizedContactFrame(2);
                dlambdaTy = dlambda * vNormalizedContactFrame(3);
                lambdaNorLenMu = mu*this.lambda(1);
                lambdaTan = [this.lambda(2) + dlambdaTx;this.lambda(3) + dlambdaTy];
                lambdaTanLen = norm(lambdaTan);
                dlambdaTan = [dlambdaTx; dlambdaTy];
				if lambdaTanLen > lambdaNorLenMu
                    dlambdaTan = lambdaTan / lambdaTanLen * lambdaNorLenMu - [this.lambda(2); this.lambda(3)];
                end
				this.lambda(2) = this.lambda(2) + dlambdaTan(1);
				this.lambda(3) = this.lambda(3) + dlambdaTan(2);
            end
            
            frictionalContactLambda = [dlambdaNor; dlambdaTan];
            dlambda = norm(frictionalContactLambda);
            if dlambda > 0
                frictionalContactNormal = [-this.nw, tx, ty] * frictionalContactLambda ./ dlambda;
                [dq1,dp1,dq2,dp2] = this.computeDx(dlambda, frictionalContactNormal);
                if this.shockProp
		            this.body1.dxJacobiShock(1:4) = this.body1.dxJacobiShock(1:4) + dq1;
                    this.body1.dxJacobiShock(5:7) = this.body1.dxJacobiShock(5:7) + dp1;
                else
		            this.body1.dxJacobi(1:4) = this.body1.dxJacobi(1:4) + dq1;
		            this.body1.dxJacobi(5:7) = this.body1.dxJacobi(5:7) + dp1;
                end
		        this.body2.dxJacobi(1:4) = this.body2.dxJacobi(1:4) + dq2;
		        this.body2.dxJacobi(5:7) = this.body2.dxJacobi(5:7) + dp2;
            end

		end

		%%
		function dlambda = solvePosDir2(this,c,nw)
			% Use the provided normal rather than normalizing
			m1 = this.body1.Mp;
			m2 = this.body2.Mp;
			I1 = this.body1.Mr;
			I2 = this.body2.Mr;
			q1 = this.body1.x(1:4);
			q2 = this.body2.x(1:4);
			nl1 = se3.qRotInv(q1,nw);
			nl2 = se3.qRotInv(q2,nw);
			rl1 = this.x1;
			rl2 = this.x2;
			rnl1 = se3.cross(rl1,nl1);
			rnl2 = se3.cross(rl2,nl2);
			w1 = (1/m1) + rnl1'*(I1.\rnl1);
			w2 = (1/m2) + rnl2'*(I2.\rnl2);
			numerator = -c;
			denominator = w1 + w2;
			dlambda = numerator/denominator;
		end

		%%
		function draw(this)
			x = this.body1.transformPoint(this.x1);
			plot3(x(1),x(2),x(3),'ro');
			x = this.body2.transformPoint(this.x2);
			plot3(x(1),x(2),x(3),'go');
		end
	end
end
