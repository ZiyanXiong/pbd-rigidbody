classdef ConCollGroundRigid < apbd.ConColl
	%ConCollGroundRigid Collision between a rigid body and the ground

	properties
		body
		xl % collision point wrt body (3x1)
		xw % collision point wrt world (3x1)
		vw % collision velocity wrt world (3x1)
		Eg % ground transformation

        contactFrame
        w1   % Generalized mass vector (3x1)
        delLinVel1 % Unit change for linear velocity matrix(3x3)
        angDelta1 % Unit change for angular velocity matrix(3x3)
        raXn  % ra X nw * sqrt(I^(-1)) matrix(3x3)

        mu
        biasCoefficient
        collision
	end

	methods
		%%
        function this = ConCollGroundRigid(body,c, collision)
			this = this@apbd.ConColl();
			this.body = body;
			this.nw = c.nw;
			this.xl = c.x1;
			this.xw = c.x2;

            this.contactFrame = zeros(3,3);
            this.w1 = zeros(3,1);
            this.raXn = zeros(3,3);
            this.delLinVel1 = zeros(3,3);
            this.angDelta1 = zeros(3,3);
            this.collision = collision;
		end

		%%
		function init(this, h, hs)
            this.d = this.body.transformPoint(this.xl) - this.xw;
            scale = min([0.8 2 * sqrt(hs / h)]);
            if this.nw' * this.d <= 0
                this.biasCoefficient = -scale / hs;
            else
                this.biasCoefficient = -1 / hs;
            end
            this.lambda = zeros(3,1);
            [tanx,tany] = apbd.ConColl.generateTangents(this.nw);
            this.contactFrame = [this.nw, tanx, tany];
            this.mu = this.body.mu;

			m1 = this.body.Mp;
			I1 = this.body.Mr;
			q1 = this.body.x0(1:4);
			rl1 = this.xl;

            for i = 1:3
                nl1 = se3.qRotInv(q1, this.contactFrame(:,i));
			    rnl1 = se3.cross(rl1,nl1);
                raXnI1 = se3.qRot(q1,(sqrt(I1).\rnl1));
			    this.w1(i) = (1/m1) +raXnI1' * raXnI1;
                this.raXn(:,i) = se3.qRot(q1,rnl1);
                
                this.delLinVel1(:,i) = this.contactFrame(:,i) / m1;
                this.angDelta1(:,i) = se3.qRot(q1,(I1.\rnl1));
            end
        end

        function layer = getLayer(this)
            layer = this.body.layer;
        end

		%%
        function solveNorPos(this, minpenetration, withSP)
            %sep = this.nw' * (this.body.transformPoint(this.xl) - this.xw0) + this.d;
            sep = this.nw' * this.body.deltaLinDt + this.raXn(:,1)' * this.body.deltaAngDt + this.nw'* this.d;
            sep = max(minpenetration,sep);
            bias = sep * this.biasCoefficient;
            %normalVel = this.nw' * this.body.computePointVel(this.xl);
            normalVel = this.nw .* this.body.v + this.body.w .* this.raXn(:,1);
            this.dlambdaNor =  bias / this.w1(1) - sum(normalVel) / this.w1(1);
            lambda = this.lambda(1) + this.dlambdaNor;
            if(lambda < 0)
                this.dlambdaNor = - this.lambda(1);
                this.collision.broken = true;
            end
            this.lambda(1) = this.lambda(1) + this.dlambdaNor;
            this.body.v = this.body.v + this.dlambdaNor * this.delLinVel1(:,1);
            this.body.w = this.body.w + this.dlambdaNor * this.angDelta1(:,1);
        end

		%%
        function solveTanVel(this, withSP)
            dlambdaTan = zeros(2,1);
            for i = 2:3
                sep = this.contactFrame(:,i)' * this.body.deltaLinDt + this.raXn(:,i)' * this.body.deltaAngDt + this.contactFrame(:,i)' * this.d;
                bias = sep * this.biasCoefficient;
                normalVel = this.contactFrame(:,i) .* this.body.v + this.body.w .* this.raXn(:,i);
                dlambdaTan(i-1) =  (bias / this.w1(i) - sum(normalVel) / this.w1(i))*0.8;
            end
            dlambdas = [0;dlambdaTan];
            %dlambdas = this.wMat \ b;
            lambdas = this.lambda + dlambdas;
            frictionRadius = this.mu * lambdas(1);
            if(norm(lambdas(2:3)) > frictionRadius)
                lambdas(2:3) = frictionRadius * lambdas(2:3) / norm(lambdas(2:3));
                dlambdas = lambdas - this.lambda; 
                this.collision.broken = true;
            end
            this.lambda = this.lambda + dlambdas;
            this.body.v = this.body.v + this.delLinVel1 * dlambdas;
            this.body.w = this.body.w + this.angDelta1 * dlambdas;
        end

        %%
        function applyLambdaSP(this)
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
