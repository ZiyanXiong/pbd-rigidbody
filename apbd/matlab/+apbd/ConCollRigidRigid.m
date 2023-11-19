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
		function update(this)
			% Assume the local positions and the normal are fixed. Compute
			% the new world positions and the new distance.
			xw1 = this.body1.transformPoint(this.x1);
			xw2 = this.body2.transformPoint(this.x2);

            q = this.body1.x0(1:4);
            p = this.body1.x0(5:7);
			xw1i = se3.qRot(q,this.x1) + p;

            q = this.body2.x0(1:4);
            p = this.body2.x0(5:7);
			xw2i = se3.qRot(q,this.x2) + p;

			% The normal stored in this object points from body1 to body2,
			% and a collision occurs if the distance is negative.
			%dval = (xw2 - xw2i) + (xw1 - xw1i);
            dval = (xw2 - xw1) - (xw2i - xw1i);
			this.d = this.nw'*dval;
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
			thresh = 0; % threshold for not fully pushing out the contact point
            
            %{
			v1w = this.body1.computePointVel(this.x1, hs);
			v2w = this.body2.computePointVel(this.x2, hs);
			v = v2w - v1w;
            this.d = hs * this.nw' * v;
            
            %this.update();

            nw = -this.nw;
            dist = (1 - thresh) * this.d;
			this.dlambdaNor = this.solvePosDir2(dist, nw);
			this.C(1) = dist;

			lambdaNor = this.lambda(1) + this.dlambdaNor;
            if lambdaNor < 0
                this.dlambdaNor = - this.lambda(1);
            end
			this.lambda(1) = this.lambda(1) + this.dlambdaNor;
            [dq1,dp1,dq2,dp2] = this.computeDx(this.dlambdaNor, nw);
		    % Save Jacobi updates
            if this.shockProp
		        this.body1.dxJacobiShock(1:4) = this.body1.dxJacobiShock(1:4) + dq1;
                this.body1.dxJacobiShock(5:7) = this.body1.dxJacobiShock(5:7) + dp1;
            else
		        this.body1.dxJacobi(1:4) = this.body1.dxJacobi(1:4) + dq1;
		        this.body1.dxJacobi(5:7) = this.body1.dxJacobi(5:7) + dp1;
            end
		    this.body2.dxJacobi(1:4) = this.body2.dxJacobi(1:4) + dq2;
		    this.body2.dxJacobi(5:7) = this.body2.dxJacobi(5:7) + dp2;
            %}

			v1w = this.body1.computePointVel(this.x1, hs);
			v2w = this.body2.computePointVel(this.x2, hs);
			v = hs * (v1w - v2w);
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
		function solveTanVel(this, hs)
			%this.update();
			mu1 = this.body1.mu;
			mu2 = this.body2.mu;
			mu = 0.5*(mu1 + mu2);
			if mu > 0 
				[tx,ty] = apbd.ConColl.generateTangents(this.nw);
                
                
				v1w = this.body1.computePointVel(this.x1, hs);
				v2w = this.body2.computePointVel(this.x2, hs);
				v = hs * (v1w - v2w);

                vt = v - this.nw * this.nw' * v;
                nt = vt ./ norm(vt);
                vt = norm(vt);

                dlambdaTan = this.solvePosDir2(vt,nt);
				this.C(2) = vt * nt' * tx;
				this.C(3) = vt * nt' * ty;
				dlambdaTx = dlambdaTan * nt' * tx;
				dlambdaTy = dlambdaTan * nt' * ty;
				% Friction limit
				lambdaNorLenMu = mu*this.lambda(1);
				lambdaTan = [this.lambda(2) + dlambdaTx;this.lambda(3) + dlambdaTy];
				lambdaTanLen = norm(lambdaTan);
                dlambdaTan = [dlambdaTx, dlambdaTy];
				if lambdaTanLen > lambdaNorLenMu
                    dlambdaTan = lambdaTan / lambdaTanLen * lambdaNorLenMu - [this.lambda(2); this.lambda(3)];
                end

                [dqTx1,dpTx1,dqTx2,dpTx2] = this.computeDx(dlambdaTan(1), tx);
                [dqTy1,dpTy1,dqTy2,dpTy2] = this.computeDx(dlambdaTan(2), ty);
				this.lambda(2) = this.lambda(2) + dlambdaTan(1);
				this.lambda(3) = this.lambda(3) + dlambdaTan(2);
				% Save Jacobi updates
                if this.shockProp
		            this.body1.dxJacobiShock(1:4) = this.body1.dxJacobiShock(1:4) + (dqTx1 + dqTy1);
                    this.body1.dxJacobiShock(5:7) = this.body1.dxJacobiShock(5:7) + (dpTx1 + dpTy1);
                else
		            this.body1.dxJacobi(1:4) = this.body1.dxJacobi(1:4) + (dqTx1 + dqTy1);
		            this.body1.dxJacobi(5:7) = this.body1.dxJacobi(5:7) + (dpTx1 + dpTy1);
                end
				this.body2.dxJacobi(1:4) = this.body2.dxJacobi(1:4) + (dqTx2 + dqTy2);
				this.body2.dxJacobi(5:7) = this.body2.dxJacobi(5:7) + (dpTx2 + dpTy2);

				%vx = hs*(v'*tx);
				%vy = hs*(v'*ty);
                
                
                %{
			    xw1 = this.body1.transformPoint(this.x1);
			    xw2 = this.body2.transformPoint(this.x2);
    
                q = this.body1.x0(1:4);
                p = this.body1.x0(5:7);
			    xw1i = se3.qRot(q,this.x1) + p;
    
                q = this.body2.x0(1:4);
                p = this.body2.x0(5:7);
			    xw2i = se3.qRot(q,this.x2) + p;
    
			    % The normal stored in this object points from body1 to body2,
			    % and a collision occurs if the distance is negative.
			    %dval = (xw2 - xw2i) + (xw1 - xw1i);
                 dval = (xw1 - xw2) - (xw1i - xw2i);
			    vx = dval'*tx;
                vy = dval'*ty;
                %}

                %{
				this.C(2) = vx;
				this.C(3) = vy;
				dlambdaTx = this.solvePosDir2(vx,tx);
				dlambdaTy = this.solvePosDir2(vy,ty);
				% Friction limit
				lambdaNorLenMu = mu*this.lambda(1);
				lambdaTan = [this.lambda(2) + dlambdaTx;this.lambda(3) + dlambdaTy];
				lambdaTanLen = norm(lambdaTan);
                dlambdaTan = [dlambdaTx, dlambdaTy];
				if lambdaTanLen > lambdaNorLenMu
                    dlambdaTan = lambdaTan / lambdaTanLen * lambdaNorLenMu - [this.lambda(2); this.lambda(3)];
                end

                [dqTx1,dpTx1,dqTx2,dpTx2] = this.computeDx(dlambdaTan(1), tx);
                [dqTy1,dpTy1,dqTy2,dpTy2] = this.computeDx(dlambdaTan(2), ty);
				this.lambda(2) = this.lambda(2) + dlambdaTan(1);
				this.lambda(3) = this.lambda(3) + dlambdaTan(2);
				% Save Jacobi updates
                if this.shockProp
		            this.body1.dxJacobiShock(1:4) = this.body1.dxJacobiShock(1:4) + (dqTx1 + dqTy1);
                    this.body1.dxJacobiShock(5:7) = this.body1.dxJacobiShock(5:7) + (dpTx1 + dpTy1);
                else
		            this.body1.dxJacobi(1:4) = this.body1.dxJacobi(1:4) + (dqTx1 + dqTy1);
		            this.body1.dxJacobi(5:7) = this.body1.dxJacobi(5:7) + (dpTx1 + dpTy1);
                end
				this.body2.dxJacobi(1:4) = this.body2.dxJacobi(1:4) + (dqTx2 + dqTy2);
				this.body2.dxJacobi(5:7) = this.body2.dxJacobi(5:7) + (dpTx2 + dpTy2);
                %}
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
