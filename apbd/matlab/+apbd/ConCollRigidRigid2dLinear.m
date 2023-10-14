classdef ConCollRigidRigid2dLinear < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
		body1
		body2
		x1 % Position wrt body 1 (3x1)
		x2 % Position wrt body 2 (3x1)
	end

	methods
		%%
        function this = ConCollRigidRigid2dLinear(body1,body2)
			this = this@apbd.ConColl();
			this.body1 = body1;
			this.body2 = body2;
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
		function update(this)
			% Assume the local positions and the normal are fixed. Compute
			% the new world positions and the new distance.
			xw1 = this.body1.transformPoint(this.x1);
			xw2 = this.body2.transformPoint(this.x2);

			[q,p] = apbd.BodyRigid2d.unproj(this.body1.x0);
			xw1i = se3.qRot(q,this.x1) + p;

			[q,p] = apbd.BodyRigid2d.unproj(this.body2.x0);
			xw2i = se3.qRot(q,this.x2) + p;

			% The normal stored in this object points from body1 to body2,
			% and a collision occurs if the distance is negative.
			%dval = (xw2 - xw2i) + (xw1 - xw1i);
            dval = (xw2 - xw1) - (xw2i - xw1i);
			this.d = this.nw'*dval;
		end

		%%
		function solveNorPos(this)
			thresh = 1e-5; % threshold for not fully pushing out the contact point
			            
            dist = (1 - thresh)* this.d;
            nw = -this.nw;

            %{
			[this.dlambdaNor,dq1,dp1,dq2,dp2] = this.solvePosDir2(dist,nw);
			this.C(1) = dist;
			this.lambda(1) = this.lambda(1) + this.dlambdaNor;
			% Save Jacobi updates
		    this.body1.dxJacobi(1:2) = this.body1.dxJacobi(1:2) + dq1(3:4);
		    this.body1.dxJacobi(3:4) = this.body1.dxJacobi(3:4) + dp1(1:2);
		    this.body2.dxJacobi(1:2) = this.body2.dxJacobi(1:2) + dq2(3:4);
		    this.body2.dxJacobi(3:4) = this.body2.dxJacobi(3:4) + dp2(1:2);
            %}

            
			%if dist <= 0.0
            if true 
            %if this.lambda(1) >= 0.0
				% Negate the normal here, since the rest of the code
				% assumes that the normal points into body 1. (Ground
				% contact uses the ground normal, which points into the
				% body.)
				this.C(1) = dist;
    			[this.dlambdaNor,dq1,dp1,dq2,dp2] = this.solvePosDir2(dist,nw);
				this.lambda(1) = this.lambda(1) + this.dlambdaNor;
				% Save Jacobi updates
			    this.body1.dxJacobi(1:2) = this.body1.dxJacobi(1:2) + dq1(3:4);
			    this.body1.dxJacobi(3:4) = this.body1.dxJacobi(3:4) + dp1(1:2);
			    this.body2.dxJacobi(1:2) = this.body2.dxJacobi(1:2) + dq2(3:4);
			    this.body2.dxJacobi(3:4) = this.body2.dxJacobi(3:4) + dp2(1:2);
			end
            
		end

		%%
		function solveTanVel(this,k,ks,hs)
			%this.update();
			mu1 = this.body1.mu;
			mu2 = this.body2.mu;
			mu = 0.5*(mu1 + mu2);
			if mu > 0 
				[tx,ty] = apbd.ConColl.generateTangents(this.nw);
				v1w = this.body1.computePointVel(this.x1,k,ks,hs);
				v2w = this.body2.computePointVel(this.x2,k,ks,hs);
				v = v1w - v2w;
				vx = hs*(v'*tx);
				vy = hs*(v'*ty);
				this.C(2) = vx;
				this.C(3) = vy;
				[dlambdaTx,dqTx1,dpTx1,dqTx2,dpTx2] = this.solvePosDir2(vx,tx);
				[dlambdaTy,dqTy1,dpTy1,dqTy2,dpTy2] = this.solvePosDir2(vy,ty);
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
                dq1 = scale*(dqTx1 + dqTy1);
                dq2 = scale*(dqTx2 + dqTy2);
                dp1 = scale*(dpTx1 + dpTy1);
                dp2 = scale*(dpTx2 + dpTy2);
			    this.body1.dxJacobi(1:2) = this.body1.dxJacobi(1:2) + dq1(3:4);
			    this.body1.dxJacobi(3:4) = this.body1.dxJacobi(3:4) + dp1(1:2);
			    this.body2.dxJacobi(1:2) = this.body2.dxJacobi(1:2) + dq2(3:4);
			    this.body2.dxJacobi(3:4) = this.body2.dxJacobi(3:4) + dp2(1:2);
			end
		end

		%%
		function [dlambda,dq1,dp1,dq2,dp2] = solvePosDir2(this,c,nw)
			% Use the provided normal rather than normalizing
			m1 = this.body1.Mp;
			m2 = this.body2.Mp;
			I1 = this.body1.Mr;
			I2 = this.body2.Mr;
			q1 = apbd.BodyRigid2d.unproj(this.body1.x);
            q2 = apbd.BodyRigid2d.unproj(this.body2.x);
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
			% Position update
			dpw = dlambda*nw;
			dp1 =  dpw/m1;
			dp2 = -dpw/m2;
			% Quaternion update
			%dpl1 = se3.qRotInv(q1,dpw);
			%dpl2 = se3.qRotInv(q2,dpw);
			%qtmp1 = [se3.qRot(q1,I1.\se3.cross(rl1,dpl1)); 0];
			%qtmp2 = [se3.qRot(q2,I2.\se3.cross(rl2,dpl2)); 0];
            qtmp1 = [I1.\se3.cross(rl1,dpw); 0];
            qtmp2 = [I2.\se3.cross(rl2,dpw); 0];
			dq1 =  0.5*se3.qMul(qtmp1,q1);
			dq2 = -0.5*se3.qMul(qtmp2,q2);
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
