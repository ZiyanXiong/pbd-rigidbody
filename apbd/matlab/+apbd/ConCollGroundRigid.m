classdef ConCollGroundRigid < apbd.ConColl
	%ConCollGroundRigid Collision between a rigid body and the ground

	properties
		body
		xl % collision point wrt body (3x1)
		xw % collision point wrt world (3x1)
		vw % collision velocity wrt world (3x1)
		Eg % ground transformation
        dt % d/h

        contactFrame
        w1   % Generalized mass vector (3x1)
        delLinVel1 % Unit change for linear velocity matrix(3x3)
        angDelta1 % Unit change for angular velocity matrix(3x3)
        raXn  % ra X nw * sqrt(I^(-1)) matrix(3x3)
        raXnI  % ra X nw * sqrt(I^(-1)) matrix(3x3)

        mu
        biasCoefficient
        collision
        dlambdas
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
            this.dlambdas = zeros(3,1);
            this.collision = collision;
		end

		%%
		function init(this, h, hs)
            this.d = this.body.transformPoint(this.xl) - this.xw;
            this.dt = this.d / h;
            this.biasCoefficient = -1 / hs;
            this.lambda = zeros(3,1);
            [tanx,tany] = apbd.ConColl.generateTangents(this.nw);
            this.contactFrame = [this.nw, tanx, tany];
            this.mu = this.body.mu;

			m1 = this.body.Mp;
			I1 = this.body.Mr;
			q1 = this.body.x(1:4);
			rl1 = this.xl;

            for i = 1:3
                nl1 = se3.qRotInv(q1, this.contactFrame(:,i));
			    rnl1 = se3.cross(rl1,nl1);
                this.raXnI(:,i) = se3.qRot(q1,(sqrt(I1).\rnl1));
			    this.w1(i) = (1/m1) + this.raXnI(:,i)' * this.raXnI(:,i);
                this.raXn(:,i) = se3.qRot(q1,rnl1);
                
                this.delLinVel1(:,i) = this.contactFrame(:,i) / m1;
                this.angDelta1(:,i) = se3.qRot(q1,(I1.\rnl1));
            end
            if isinf(m1)
                this.w1 = ones(3,1);
            end
        end

        function layer = getLayer(this)
            layer = this.body.layer;
        end

        %%
        function Cs = evalCs(this)
            Cs = this.contactFrame' * this.body.v + this.raXn' * this.body.w + this.contactFrame'* this.dt;
        end

        %%
        function applyLambda(this, dlambdas)
            this.lambda = this.lambda + dlambdas;
            this.body.v = this.body.v + this.delLinVel1 * dlambdas;
            this.body.w = this.body.w + this.angDelta1 * dlambdas;
        end

		%%
        function solveNorPos(this, withSP)
            if(~withSP)
                sep = this.nw' * this.body.deltaLinDt + this.raXn(:,1)' * this.body.deltaAngDt + this.nw'* this.d;
                bias = sep * this.biasCoefficient;
                %normalVel = this.nw' * this.body.computePointVel(this.xl);
                normalVel = this.nw .* this.body.v + this.body.w .* this.raXn(:,1);
                this.dlambdas(1) =  bias / this.w1(1) - sum(normalVel) / this.w1(1);
                lambda = this.lambda(1) + this.dlambdas(1);
                if(lambda < 0)
                    this.dlambdas(1) = - this.lambda(1);
                    %this.collision.broken = true;
                end
                this.lambda(1) = this.lambda(1) + this.dlambdas(1);
                this.body.v = this.body.v + this.dlambdas(1) * this.delLinVel1(:,1);
                this.body.w = this.body.w + this.dlambdas(1) * this.angDelta1(:,1);
            else
                sep = this.nw' * this.body.deltaLinDt + this.raXn(:,1)' * this.body.deltaAngDt + this.nw'* this.dt;
                bias = -sep;
                %normalVel = this.nw' * this.body.computePointVel(this.xl);
                normalVel = this.nw .* this.body.v + this.body.w .* this.raXn(:,1);
                this.dlambdas(1) =  bias / this.w1(1) - sum(normalVel) / this.w1(1);
                lambda = this.lambda(1) + this.dlambdas(1);
                if(lambda < 0)
                    this.dlambdas(1) = - this.lambda(1);
                    this.collision.broken = true;
                end
                this.lambda(1) = this.lambda(1) + this.dlambdas(1);
                this.body.v = this.body.v + this.dlambdas(1) * this.delLinVel1(:,1);
                this.body.w = this.body.w + this.dlambdas(1) * this.angDelta1(:,1);
            end
        end

		%%
        function solveTanPos(this, withSP)
            dlambdaTan = zeros(2,1);
            if(~withSP)
                for i = 2:3
                    sep = this.contactFrame(:,i)' * this.body.deltaLinDt + this.raXn(:,i)' * this.body.deltaAngDt + this.contactFrame(:,i)' * this.d;
                    bias = sep * this.biasCoefficient;
                    normalVel = this.contactFrame(:,i) .* this.body.v + this.body.w .* this.raXn(:,i);
                    dlambdaTan(i-1) =  (bias / this.w1(i) - sum(normalVel) / this.w1(i));
                end
                dlambdaTan = [0;dlambdaTan];
                %dlambdas = this.wMat \ b;
                lambdas = this.lambda + dlambdaTan;
                frictionRadius = this.mu * lambdas(1);
                if(norm(lambdas(2:3)) > frictionRadius)
                    lambdas(2:3) = frictionRadius * lambdas(2:3) / norm(lambdas(2:3));
                    dlambdaTan = lambdas - this.lambda; 
                    %this.collision.broken = true;
                end
                this.lambda = this.lambda + dlambdaTan;
                this.body.v = this.body.v + this.delLinVel1 * dlambdaTan;
                this.body.w = this.body.w + this.angDelta1 * dlambdaTan;
                this.dlambdas(2:3) = dlambdaTan(2:3);
            else
                for i = 2:3
                    sep = this.contactFrame(:,i)' * this.body.deltaLinDt + this.raXn(:,i)' * this.body.deltaAngDt + this.contactFrame(:,i)' * this.dt;
                    bias = -sep;
                    normalVel = this.contactFrame(:,i) .* this.body.v + this.body.w .* this.raXn(:,i);
                    dlambdaTan(i-1) =  (bias / this.w1(i) - sum(normalVel) / this.w1(i));
                end
                dlambdaTan = [0;dlambdaTan];
                %dlambdas = this.wMat \ b;
                lambdas = this.lambda + dlambdaTan;
                frictionRadius = this.mu * lambdas(1);
                if(norm(lambdas(2:3)) > frictionRadius)
                    lambdas(2:3) = frictionRadius * lambdas(2:3) / norm(lambdas(2:3));
                    dlambdaTan = lambdas - this.lambda; 
                    this.collision.broken = true;
                end
                this.lambda = this.lambda + dlambdaTan;
                this.body.v = this.body.v + this.delLinVel1 * dlambdaTan;
                this.body.w = this.body.w + this.angDelta1 * dlambdaTan;
                this.dlambdas(2:3) = dlambdaTan(2:3);
            end
        end

		%%
        function solveNorVel(this, substeps)
            sep = this.nw' * this.body.deltaLinDt + this.raXn(:,1)' * this.body.deltaAngDt;
            bias = -sep * substeps;
            %normalVel = this.nw' * this.body.computePointVel(this.xl);
            normalVel = this.nw .* this.body.v + this.body.w .* this.raXn(:,1);
            this.dlambdas(1) =  bias / this.w1(1) - sum(normalVel) / this.w1(1);
            lambda = this.lambda(1) + this.dlambdas(1);
            if(lambda < 0)
                this.dlambdas(1) = - this.lambda(1);
                %this.collision.broken = true;
            end
            this.lambda(1) = this.lambda(1) + this.dlambdas(1);
            this.body.v = this.body.v + this.dlambdas(1) * this.delLinVel1(:,1);
            this.body.w = this.body.w + this.dlambdas(1) * this.angDelta1(:,1);
        end

		%%
        function solveTanVel(this,substeps)
            dlambdaTan = zeros(2,1);
            for i = 2:3
                sep = this.contactFrame(:,i)' * this.body.deltaLinDt + this.raXn(:,i)' * this.body.deltaAngDt;
                bias = -sep * substeps;
                normalVel = this.contactFrame(:,i) .* this.body.v + this.body.w .* this.raXn(:,i);
                dlambdaTan(i-1) =  (bias / this.w1(i) - sum(normalVel) / this.w1(i));
            end
            dlambdaTan = [0;dlambdaTan];
            %dlambdas = this.wMat \ b;
            lambdas = this.lambda + dlambdaTan;
            frictionRadius = this.mu * lambdas(1);
            if(norm(lambdas(2:3)) > frictionRadius)
                lambdas(2:3) = frictionRadius * lambdas(2:3) / norm(lambdas(2:3));
                dlambdaTan = lambdas - this.lambda; 
                %this.collision.broken = true;
            end
            this.lambda = this.lambda + dlambdaTan;
            this.body.v = this.body.v + this.delLinVel1 * dlambdaTan;
            this.body.w = this.body.w + this.angDelta1 * dlambdaTan;
            this.dlambdas(2:3) = dlambdaTan(2:3);
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
