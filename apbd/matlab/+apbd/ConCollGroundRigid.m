classdef ConCollGroundRigid < apbd.ConColl
	%ConCollGroundRigid Collision between a rigid body and the ground

	properties
		body
		xl % collision point wrt body (3x1)
		xw % collision point wrt world (3x1)
		vw % collision velocity wrt world (3x1)
		Eg % ground transformation
	end

	methods
		%%
		function this = ConCollGroundRigid(body,Eg)
			this = this@apbd.ConColl();
			this.body = body;
            this.bodies = {body};
			this.Eg = Eg;
			this.xl = zeros(3,1);
			this.xw = zeros(3,1);
			this.vw = zeros(3,1);
            this.ground = true;
		end

		%%
		function init(this) %#ok<MANU>
			% Do nothing
        end

        %%
        function applyJacobi(this)
            this.body.applyJacobi();
        end
        
		%%
		function update(this)
            
			this.xw = this.body.transformPoint(this.xl);
			xg = this.Eg\[this.xw;1];

            q = this.body.x0(1:4);
            p = this.body.x0(5:7);
			xwi = se3.qRot(q,this.xl) + p;
            xgi = this.Eg\[xwi;1];
			this.d = xg(3) - xgi(3);
        end

        %%
        function [dq,dp] = computeDx(this, dlambda, nw)
			m1 = this.body.Mp;
			I1 = this.body.Mr;
			% Position update
			dpw = dlambda*nw;
			dp = dpw/m1;
			% Quaternion update
            q1 = this.body.x1_0(1:4);
			dpl1 = se3.qRotInv(q1,dpw);
			qtmp1 = [se3.qRot(q1, I1.\se3.cross(this.xl,dpl1)); 0];
            %qtmp1 = [I1.\se3.cross(rl1,dpl1); 0];
			%dq = se3.qMul(sin(0.5*qtmp1),q1);
            dq = 0.5 * se3.qMul(qtmp1,q1);
		end

		%%
		function solveNorPos(this, hs)
			thresh = 0; % threshold for not fully pushing out the contact point
			%dist = (1 - thresh)*this.d;
			%fprintf('%d ',dist < 0);
            
            %{
            v = this.body.computePointVel(this.xl,hs);
            this.d = hs * this.nw' * v;
            
            
            %this.update();

			dist = (1 - thresh)*this.d;
            this.dlambdaNor = this.solvePosDir1(dist,this.nw);
			this.C(1) = dist;
			lambdaNor = this.lambda(1) + this.dlambdaNor;
            if lambdaNor < 0
                this.dlambdaNor = - this.lambda(1);
            end
            this.lambda(1) = this.lambda(1) + this.dlambdaNor;
            [dq,dp] = this.computeDx(this.dlambdaNor, this.nw);
			% Save Jacobi updates
			this.body.dxJacobi(1:4) = this.body.dxJacobi(1:4) + dq;
			this.body.dxJacobi(5:7) = this.body.dxJacobi(5:7) + dp;
            %}

            v = hs * this.body.computePointVel(this.xl, hs);
            vNorm = norm(v);
            vNormalized = v ./ vNorm;
			tx = this.Eg(1:3,1);
			ty = this.Eg(1:3,2);
            vNormalizedContactFrame = [this.nw'; tx' ; ty'] * vNormalized;

            dlambda = this.solvePosDir1(norm(v), vNormalized);
            this.C = vNorm * vNormalizedContactFrame;

            dlambdaNor = dlambda * vNormalizedContactFrame(1);
            lambdaNor = this.lambda(1) + dlambdaNor;
            if lambdaNor < 0
                dlambdaNor  = - this.lambda(1);
            end
            this.lambda(1) = this.lambda(1) + dlambdaNor;
            mu = this.body.mu;
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
                frictionalContactNormal = [this.nw, tx, ty] * frictionalContactLambda ./ dlambda;
                [dq,dp] = this.computeDx(dlambda, frictionalContactNormal);
			    this.body.dxJacobi(1:4) = this.body.dxJacobi(1:4) + dq;
			    this.body.dxJacobi(5:7) = this.body.dxJacobi(5:7) + dp;
            end
		end

		%%
		function solveTanVel(this, hs)
			mu = this.body.mu;
			if mu > 0 
				%[tx,ty] = apbd.ConColl.generateTangents(this.nw);
				tx = this.Eg(1:3,1);
				ty = this.Eg(1:3,2);
				v = hs * this.body.computePointVel(this.xl, hs);

                vt = v - this.nw * this.nw' * v;
                nt = vt ./ norm(vt);
                vt = norm(vt);
                
                dlambdaTan = this.solvePosDir1(vt,nt);
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
                dlambdaNt = (dlambdaTan(1) * tx + dlambdaTan(2) * ty)' * nt;
                %[dqTx,dpTx] = this.computeDx(dlambdaTan(1), tx);
                %[dqTy,dpTy] = this.computeDx(dlambdaTan(2), ty);
                [dqTan,dpTan] = this.computeDx(dlambdaNt, nt);
				this.lambda(2) = this.lambda(2) + dlambdaTan(1);
				this.lambda(3) = this.lambda(3) + dlambdaTan(2);
				% Save Jacobi updates
				%this.body.dxJacobi(1:4) = this.body.dxJacobi(1:4) + (dqTx + dqTy);
				%this.body.dxJacobi(5:7) = this.body.dxJacobi(5:7) + (dpTx + dpTy);
				this.body.dxJacobi(1:4) = this.body.dxJacobi(1:4) + dqTan;
				this.body.dxJacobi(5:7) = this.body.dxJacobi(5:7) + dpTan;
                

				%vx = hs*(v'*tx);
				%vy = hs*(v'*ty);
                
                
                %{
			    this.xw = this.body.transformPoint(this.xl);
                q = this.body.x0(1:4);
                p = this.body.x0(5:7);
			    xwi = se3.qRot(q,this.xl) + p;
                dval = this.xw - xwi;
			    vx = dval'*tx;
                vy = dval'*ty;
                %}
                
                %{
				this.C(2) = vx;
				this.C(3) = vy;
				dlambdaTx = this.solvePosDir1(vx,tx);
				dlambdaTy = this.solvePosDir1(vy,ty);
				% Friction limit
				lambdaNorLenMu = mu*this.lambda(1);
				lambdaTan = [this.lambda(2) + dlambdaTx;this.lambda(3) + dlambdaTy];
				lambdaTanLen = norm(lambdaTan);
                dlambdaTan = [dlambdaTx, dlambdaTy];
				if lambdaTanLen > lambdaNorLenMu
                    dlambdaTan = lambdaTan / lambdaTanLen * lambdaNorLenMu - [this.lambda(2); this.lambda(3)];
                end

                [dqTx,dpTx] = this.computeDx(dlambdaTan(1), tx);
                [dqTy,dpTy] = this.computeDx(dlambdaTan(2), ty);
				this.lambda(2) = this.lambda(2) + dlambdaTan(1);
				this.lambda(3) = this.lambda(3) + dlambdaTan(2);
				% Save Jacobi updates
				this.body.dxJacobi(1:4) = this.body.dxJacobi(1:4) + (dqTx + dqTy);
				this.body.dxJacobi(5:7) = this.body.dxJacobi(5:7) + (dpTx + dpTy);
                %}
			end
		end

		%%
		function dlambda = solvePosDir1(this,c,nw)
			% Use the provided normal rather than normalizing
			m1 = this.body.Mp;
			I1 = this.body.Mr;
			q1 = this.body.x(1:4);
			nl1 = se3.qRotInv(q1,nw);
			rl1 = this.xl;
			rnl1 = se3.cross(rl1,nl1);
			w1 = (1/m1) + rnl1'*(I1.\rnl1);
			numerator = -c;
			denominator = w1;
			dlambda = numerator/denominator;
		end

		%%
		function draw(this)
			x = this.body.transformPoint(this.xl);
			plot3(x(1),x(2),x(3),'go');
			x = this.xw;
			plot3(x(1),x(2),x(3),'ro');
			x = [this.xw(1:3),this.xw(1:3)+this.s*this.nw(1:3)];
			plot3(x(1,:),x(2,:),x(3,:),'r-');
		end
	end
end
