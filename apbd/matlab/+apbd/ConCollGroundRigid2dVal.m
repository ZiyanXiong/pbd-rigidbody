classdef ConCollGroundRigid2dVal < apbd.ConColl
	%ConCollGroundRigid2d Collision between a rigid body and the ground

	properties
		body
		xl % collision point wrt body (3x1)
		xw % collision point wrt world (3x1)
		vw % collision velocity wrt world (3x1)
		Eg % ground transformation
	end

	methods
		%%
		function this = ConCollGroundRigid2dVal(body,Eg)
			this = this@apbd.ConColl();
			this.body = body;
			this.Eg = Eg;
			this.xl = zeros(3,1);
			this.xw = zeros(3,1);
			this.vw = zeros(3,1);
		end

		%%
		function init(this) %#ok<MANU>
			% Do nothing
		end

		%%
		function update(this)
			this.xw = this.body.transformPoint(this.xl);
			xg = this.Eg\[this.xw;1];

			[q,p] = apbd.BodyRigid2d.unproj(this.body.x0);
			xwi = se3.qRot(q,this.xl) + p;
            xgi = this.Eg\[xwi;1];
			this.d = xg(3) - xgi(3);
		end

		%%
		function solveNorPos(this)
			thresh = 1e-5; % threshold for not fully pushing out the contact point

            %{
            dist = (1 - thresh) * this.body.transformPoint(this.xl);
            dist = dist(2);
			this.C(1) = dist;
			[this.dlambdaNor,dq,dp] = this.solvePosDir1(dist,this.nw);
			this.lambda(1) = this.lambda(1) + this.dlambdaNor;
			% Save Jacobi updates
			this.body.dxJacobi(1:2) = this.body.dxJacobi(1:2) + dq(3:4);
			this.body.dxJacobi(3:4) = this.body.dxJacobi(3:4) + dp(1:2);
            %}
            
            
			dist = (1 - thresh)*this.d;

            %{
            [this.dlambdaNor,dq,dp] = this.solvePosDir1(dist,this.nw);

			this.C(1) = dist;
			this.lambda(1) = this.lambda(1) + this.dlambdaNor;
			% Save Jacobi updates
			this.body.dxJacobi(1:2) = this.body.dxJacobi(1:2) + dq(3:4);
			this.body.dxJacobi(3:4) = this.body.dxJacobi(3:4) + dp(1:2);
            %}

            %if dist <= 0.0
            if true
			%if this.lambda(1) >= 0.0
				this.C(1) = dist;
                [this.dlambdaNor,dq,dp] = this.solvePosDir1(dist,this.nw);
				this.lambda(1) = this.lambda(1) + this.dlambdaNor;
				% Save Jacobi updates
				this.body.dxJacobi(1:2) = this.body.dxJacobi(1:2) + dq(3:4);
				this.body.dxJacobi(3:4) = this.body.dxJacobi(3:4) + dp(1:2);
			end
            
		end

		%%
		function solveTanVel(this,k,ks,hs)
			mu = this.body.mu;
			if mu > 0 
				%[tx,ty] = apbd.ConColl.generateTangents(this.nw);
				tx = this.Eg(1:3,1);
				ty = this.Eg(1:3,2);
				v = this.body.computePointVel(this.xl,k,ks,hs);
				vx = hs*(v'*tx);
				vy = hs*(v'*ty);
				this.C(2) = vx;
				this.C(3) = vy;
				[dlambdaTx,dqTx,dpTx] = this.solvePosDir1(vx,tx);
				[dlambdaTy,dqTy,dpTy] = this.solvePosDir1(vy,ty);
				% Friction limit
				dlambdaNorLenMu = mu*this.dlambdaNor;
				dlambdaTan = [dlambdaTx;dlambdaTy];
				dlambdaTanLen = norm(dlambdaTan);
				scale = 1;
				if dlambdaTanLen > dlambdaNorLenMu
					scale = dlambdaNorLenMu/dlambdaTanLen;
				end
				this.lambda(2) = this.lambda(2) + scale*dlambdaTx;
				this.lambda(3) = this.lambda(3) + scale*dlambdaTy;
				% Save Jacobi updates
				this.body.dxJacobi(1:2) = this.body.dxJacobi(1:2) + scale*(dqTx(3:4) + dqTy(3:4));
				this.body.dxJacobi(3:4) = this.body.dxJacobi(3:4) + scale*(dpTx(1:2) + dpTy(1:2));
			end
		end

		%%
		function [dlambda,dq1,dp1] = solvePosDir1(this,c,nw)
			% Use the provided normal rather than normalizing
			m1 = this.body.Mp;
			I1 = this.body.Mr;
			q1 = apbd.BodyRigid2d.unproj(this.body.x);
			nl1 = se3.qRotInv(q1,nw);
			rl1 = this.xl;
			rnl1 = se3.cross(rl1,nl1);
			w1 = (1/m1) + rnl1'*(I1.\rnl1);
			numerator = -c;
			denominator = w1;
			dlambda = numerator/denominator;
			% Position update
			dpw = dlambda*nw;
			dp1 = dpw/m1;
			% Quaternion update
			dpl1 = se3.qRotInv(q1,dpw);
			qtmp1 = [se3.qRot(q1,I1.\se3.cross(rl1,dpl1)); 0];
			dq1 = 0.5*se3.qMul(qtmp1,q1);
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
