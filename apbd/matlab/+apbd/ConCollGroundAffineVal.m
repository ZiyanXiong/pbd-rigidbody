classdef ConCollGroundAffineVal < apbd.ConColl
	%ConCollGroundAffine Collision between an affine body and the ground

	properties
		affine
		xl % collision point wrt body (3x1)
		xw % collision point wrt world (3x1)
		vw % collision velocity wrt world (3x1)
        nl % collision normal wrt to x0 (3x1)
		Eg % ground transformation
	end

	methods
		%%
		function this = ConCollGroundAffineVal(affine,Eg)
			this = this@apbd.ConColl();
			this.affine = affine;
			this.Eg = Eg;
			this.xl = zeros(3,1);
			this.xw = zeros(3,1);
			this.vw = zeros(3,1);
            this.nl = zeros(3,1);
		end

		%%
		function init(this) %#ok<MANU>
			% Do nothing
		end

		%%
		function update(this)
			this.xw = this.affine.transformPoint(this.xl);
			xg = this.Eg\[this.xw;1];

			xwi = zeros(3,1);
			xwi(1) = this.affine.x0(1:3)'*this.xl + this.affine.x0(10);
			xwi(2) = this.affine.x0(4:6)'*this.xl + this.affine.x0(11);
			xwi(3) = this.affine.x0(7:9)'*this.xl + this.affine.x0(12);
            xgi = this.Eg\[xwi;1];
			this.d = xg(3) - xgi(3);
		end

		%%
		function solveNorPos(this)
			thresh = 1e-5; % threshold for not fully pushing out the contact point
			dist = (1 - thresh)*this.d;
			%if dist < 0
            if true
				v = this.Eg(1:3,3);
                this.C(1) = dist;
				[dlambda,dx] = this.solvePos(dist, v);
				this.dlambdaNor = dlambda;
				this.lambda(1) = this.lambda(1) + dlambda;
				% Save Jacobi updates
				this.affine.dxJacobi = this.affine.dxJacobi + dx;
			end
		end

		%%
		function solveTanVel(this,k,ks,hs)
			mu = this.affine.mu;
			if mu > 0 && this.dlambdaNor > 0
				Rg = this.Eg(1:3,1:3);
				v = this.affine.computePointVel(this.xl,k,ks,hs);
				tx = Rg(:,1);
				ty = Rg(:,2);
				vx = hs*(v'*tx);
				vy = hs*(v'*ty);
				this.C(2) = vx;
				this.C(3) = vy;
				[dlambdaTx,dxTx] = this.solvePos(vx,tx);
				[dlambdaTy,dxTy] = this.solvePos(vy,ty);
				% Friction limit
				dlambdaNorLenMu = mu*this.dlambdaNor;
				dlambdaTan = [dlambdaTx;dlambdaTy];
				dlambdaTanLen = norm(dlambdaTan);
				scale = 1;
				if dlambdaTanLen > dlambdaNorLenMu
					scale = dlambdaNorLenMu/dlambdaTanLen;
				end
				dx = scale*(dxTx + dxTy);
				dlambdaTx = scale*dlambdaTx;
				dlambdaTy = scale*dlambdaTy;
				this.lambda(2) = this.lambda(2) + dlambdaTx;
				this.lambda(3) = this.lambda(3) + dlambdaTy;
				% Save Jacobi updates
				this.affine.dxJacobi = this.affine.dxJacobi + dx;
			end
		end

		%%
		function [dlambda,dx] = solvePos(this,c,v)
			numerator = -c;
			Wa = this.affine.Wa;
			Wp = this.affine.Wp;
			%J = apbd.BodyAffine.jacobian(this.xl);
			%dC = v'*J;
			%dC = [v(1)*this.xl; v(2)*this.xl; v(3)*this.xl; v]'; % row vector
			%WdC = [Wa.*dC(1:3)'; Wa.*dC(4:6)'; Wa.*dC(7:9)'; Wp*dC(10:12)']; % column vector
			dCx = v(1)*this.xl;
			dCy = v(2)*this.xl;
			dCz = v(3)*this.xl;
			dCp = v;
			WdCx = Wa.*dCx;
			WdCy = Wa.*dCy;
			WdCz = Wa.*dCz;
			WdCp = Wp*v;
			denominator = dCx'*WdCx + dCy'*WdCy + dCz'*WdCz + dCp'*WdCp;
			%denominator = dC*WdC;
			dlambda = numerator/denominator;
			dx = dlambda*[WdCx; WdCy; WdCz; WdCp];
		end

		%%
		function draw(this)
			x = this.affine.transformPoint(this.xl);
			plot3(x(1),x(2),x(3),'go');
			x = this.xw;
			plot3(x(1),x(2),x(3),'ro');
			x = [this.xw(1:3),this.xw(1:3)+this.s*this.nw(1:3)];
			plot3(x(1,:),x(2,:),x(3,:),'r-');
		end
	end
end
