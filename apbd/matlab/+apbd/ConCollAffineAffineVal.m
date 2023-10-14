classdef ConCollAffineAffineVal < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
		affine1
		affine2
		x1 % Position wrt body 1 (3x1)
		x2 % Position wrt body 2 (3x1)
	end

	methods
		%%
        function this = ConCollAffineAffineVal(affine1,affine2)
			this = this@apbd.ConColl();
			this.affine1 = affine1;
			this.affine2 = affine2;
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
			xw1 = this.affine1.transformPoint(this.x1);
			xw2 = this.affine2.transformPoint(this.x2);


			xw1i = zeros(3,1);
			xw1i(1) = this.affine1.x0(1:3)'*this.x1 + this.affine1.x0(10);
			xw1i(2) = this.affine1.x0(4:6)'*this.x1 + this.affine1.x0(11);
			xw1i(3) = this.affine1.x0(7:9)'*this.x1 + this.affine1.x0(12);

			xw2i = zeros(3,1);
			xw2i(1) = this.affine2.x0(1:3)'*this.x2 + this.affine2.x0(10);
			xw2i(2) = this.affine2.x0(4:6)'*this.x2 + this.affine2.x0(11);
			xw2i(3) = this.affine2.x0(7:9)'*this.x2 + this.affine2.x0(12);

			% The normal stored in this object points from body1 to body2,
			% and a collision occurs if the distance is negative.
			%dval = (xw2 - xw2i) + (xw1 - xw1i);
            dval = (xw2 - xw1) - (xw2i - xw1i);
			this.d = this.nw'*dval;
		end

		%%
		function solveNorPos(this)
			thresh = 1e-5; % threshold for not fully pushing out the contact point
			
            dist = (1 - thresh) * this.d;
            nw = -this.nw;

			%if dist < 0
            if true
				this.C(1) = dist;
				[dlambda,dx1, dx2] = this.solvePos2(dist,nw);
				this.dlambdaNor = dlambda;
				this.lambda(1) = this.lambda(1) + dlambda;
				% Save Jacobi updates
				this.affine1.dxJacobi = this.affine1.dxJacobi + dx1;
                this.affine2.dxJacobi = this.affine2.dxJacobi + dx2;
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
				this.body1.dxJacobi(1:4) = this.body1.dxJacobi(1:4) + scale*(dqTx1 + dqTy1);
				this.body2.dxJacobi(1:4) = this.body2.dxJacobi(1:4) + scale*(dqTx2 + dqTy2);
				this.body1.dxJacobi(5:7) = this.body1.dxJacobi(5:7) + scale*(dpTx1 + dpTy1);
				this.body2.dxJacobi(5:7) = this.body2.dxJacobi(5:7) + scale*(dpTx2 + dpTy2);
			end
		end

		%%
        function [dlambda,dx1,dx2] = solvePos2(this,c,v)
			numerator = -c;
			Wa = this.affine1.Wa;
			Wp = this.affine1.Wp;
			dCx = v(1)*this.x1;
			dCy = v(2)*this.x1;
			dCz = v(3)*this.x1;
			dCp = v;
			WdCx = Wa.*dCx;
			WdCy = Wa.*dCy;
			WdCz = Wa.*dCz;
			WdCp = Wp*v;
			w1 = dCx'*WdCx + dCy'*WdCy + dCz'*WdCz + dCp'*WdCp;
            dx1 = [WdCx; WdCy; WdCz; WdCp];

			Wa = this.affine2.Wa;
			Wp = this.affine2.Wp;
			dCx = v(1)*this.x2;
			dCy = v(2)*this.x2;
			dCz = v(3)*this.x2;
			dCp = v;
			WdCx = Wa.*dCx;
			WdCy = Wa.*dCy;
			WdCz = Wa.*dCz;
			WdCp = Wp*v;
			w2 = dCx'*WdCx + dCy'*WdCy + dCz'*WdCz + dCp'*WdCp;
            dx2 = [WdCx; WdCy; WdCz; WdCp];

			denominator = w1 + w2;
			dlambda = numerator/denominator;

			dx1 = dlambda*dx1;
            dx2 = -dlambda*dx2;
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
