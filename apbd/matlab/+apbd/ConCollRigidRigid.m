classdef ConCollRigidRigid < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
		body1
		body2
		x1 % Position wrt body 1 (3x1)
		x2 % Position wrt body 2 (3x1)

        contactFrame
        w1   % Generalized mass vector (3x1)
        delLinVel1 % Unit change for linear velocity matrix(3x3)
        angDelta1 % Unit change for angular velocity matrix(3x3)
        raXn1  % ra X nw * sqrt(I^(-1)) matrix(3x3)

        w2   % Generalized mass vector (3x1)
        delLinVel2 % Unit change for linear velocity matrix(3x3)
        angDelta2 % Unit change for angular velocity matrix(3x3)
        raXn2  % ra X nw * sqrt(I^(-1)) matrix(3x3)

        mu
        biasCoefficient
        collision
	end

	methods
		%%
        function this = ConCollRigidRigid(body1,body2, c, collision)
			this = this@apbd.ConColl();
			this.body1 = body1;
			this.body2 = body2;
			this.nw = -c.nw;
			this.x1 = c.x1;
			this.x2 = c.x2;

            this.contactFrame = zeros(3,3);

            this.w1 = zeros(3,1);
            this.raXn1 = zeros(3,3);
            this.delLinVel1 = zeros(3,3);
            this.angDelta1 = zeros(3,3);

            this.w2 = zeros(3,1);
            this.raXn2 = zeros(3,3);
            this.delLinVel2 = zeros(3,3);
            this.angDelta2 = zeros(3,3);

            this.collision = collision;
		end

		%%
		function init(this, h, hs) 
            % Contact normal always point form body2 to body1
            if(this.body1.layer < this.body2.layer)
                temp = this.body1;
                this.body1 = this.body2;
                this.body2 = temp;
                temp = this.x1;
                this.x1 = this.x2;
                this.x2 = temp;
                this.nw = - this.nw;
            end
            % Now we can assume body1 is always above body2

            this.d = this.body1.transformPoint(this.x1) - this.body2.transformPoint(this.x2);
            scale = min([0.8 2 * sqrt(hs / h)]);
            if this.nw' * this.d <= 0
                this.biasCoefficient = -scale / hs;
            else
                this.biasCoefficient = -1 / hs;
            end

            this.lambda = zeros(3,1);
            [tanx,tany] = apbd.ConColl.generateTangents(this.nw);
            this.contactFrame = [this.nw, tanx, tany];
            this.mu = 0.5 * (this.body1.mu + this.body2.mu);

			m1 = this.body1.Mp;
			I1 = this.body1.Mr;
			q1 = this.body1.x0(1:4);
			rl1 = this.x1;

			m2 = this.body2.Mp;
			I2 = this.body2.Mr;
			q2 = this.body2.x0(1:4);
			rl2 = this.x2;

            for i = 1:3
                nl1 = se3.qRotInv(q1, this.contactFrame(:,i));
			    rnl1 = se3.cross(rl1,nl1);
                raXnI1 = se3.qRot(q1,(sqrt(I1).\rnl1));
                this.raXn1(:,i) = se3.qRot(q1,rnl1);
			    this.w1(i) = (1/m1) + raXnI1' * raXnI1;
                this.delLinVel1(:,i) = this.contactFrame(:,i) / m1;
                this.angDelta1(:,i) = se3.qRot(q1,(I1.\rnl1));
    
                nl2 = se3.qRotInv(q2, this.contactFrame(:,i));
                rnl2 = se3.cross(rl2,nl2);
                raXnI2 = se3.qRot(q2,(sqrt(I2).\rnl2));
                this.raXn2(:,i) = se3.qRot(q2,rnl2);
			    this.w2(i) = (1/m2) + raXnI2' * raXnI2;
                this.delLinVel2(:,i) = this.contactFrame(:,i) / m2;
                this.angDelta2(:,i) = se3.qRot(q2,(I2.\rnl2));
            end
        end

        %%
        function layer = getLayer(this)
            layer = this.body1.layer + this.body2.layer;
        end

		%%
		function solveNorPos(this, minpenetration, withSP)
            if(~withSP)
                sep = this.nw' * this.body1.deltaLinDt + this.raXn1(:,1)' * this.body1.deltaAngDt + this.nw' * this.d;
                sep = sep - (this.nw' * this.body2.deltaLinDt + this.raXn2(:,1)' * this.body2.deltaAngDt);
                sep = max(minpenetration,sep);
                bias = sep * this.biasCoefficient;
                %normalVel = this.nw' * this.body.computePointVel(this.xl);
                normalVel = this.nw .* this.body1.v + this.body1.w .* this.raXn1(:,1);
                normalVel = normalVel - (this.nw .* this.body2.v + this.body2.w .* this.raXn2(:,1));
                this.dlambdaNor =  bias / (this.w1(1) + this.w2(1)) - sum(normalVel) / (this.w1(1) + this.w2(1));
                lambda = this.lambda(1) + this.dlambdaNor;
                if(lambda < 0)
                    this.dlambdaNor = - this.lambda(1);
                    this.collision.broken = true;
                end
                this.lambda(1) = this.lambda(1) + this.dlambdaNor;
                this.body1.v = this.body1.v + this.dlambdaNor * this.delLinVel1(:,1);
                this.body1.w = this.body1.w + this.dlambdaNor * this.angDelta1(:,1);
                this.body2.v = this.body2.v - this.dlambdaNor * this.delLinVel2(:,1);
                this.body2.w = this.body2.w - this.dlambdaNor * this.angDelta2(:,1);
            else
                sep = this.nw' * this.body1.deltaLinDt + this.raXn1(:,1)' * this.body1.deltaAngDt + this.nw' * this.d;
                sep = max(minpenetration,sep);
                bias = sep * this.biasCoefficient;
                %normalVel = this.nw' * this.body.computePointVel(this.xl);
                
                normalVel = this.nw .* this.body1.v + this.body1.w .* this.raXn1(:,1);
                this.dlambdaNor =  bias / (this.w1(1)) - sum(normalVel) / (this.w1(1));
                %this.dlambdaNor = - sum(normalVel) / (this.w1(1));
                lambda = this.lambda(1) + this.dlambdaNor;
                if(lambda < 0)
                    this.dlambdaNor = - this.lambda(1);
                    %this.collision.broken = true;
                end
                this.lambda(1) = this.lambda(1) + this.dlambdaNor;
                this.body1.v = this.body1.v + this.dlambdaNor * this.delLinVel1(:,1);
                this.body1.w = this.body1.w + this.dlambdaNor * this.angDelta1(:,1);
                this.dlambdaSP(1) = this.dlambdaSP(1) + this.dlambdaNor;
            end
        end

		%%
        function solveTanVel(this, withSP)
            dlambdaTan = zeros(2,1);
            if(~withSP)
                for i = 2:3
                    %sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXnI1(:,i)' * this.body1.deltaAngDt;
                    sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXn1(:,i)' * this.body1.deltaAngDt + this.contactFrame(:,i)' * this.d;
                    sep = sep - (this.contactFrame(:,i)' * this.body2.deltaLinDt + this.raXn2(:,i)' * this.body2.deltaAngDt);
                    bias = sep * this.biasCoefficient;
                    normalVel = this.contactFrame(:,i) .* this.body1.v + this.body1.w .* this.raXn1(:,i);
                    normalVel = normalVel - (this.contactFrame(:,i) .* this.body2.v + this.body2.w .* this.raXn2(:,i));
                    dlambdaTan(i-1) =  (bias / (this.w1(i) + this.w2(i)) - sum(normalVel) / (this.w1(i) + this.w2(i)))*0.8;
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
                this.body1.v = this.body1.v +  this.delLinVel1 * dlambdas;
                this.body1.w = this.body1.w +  this.angDelta1 * dlambdas;
                this.body2.v = this.body2.v - this.delLinVel2 * dlambdas;
                this.body2.w = this.body2.w - this.angDelta2 * dlambdas;
            else
                for i = 2:3
                    %sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXnI1(:,i)' * this.body1.deltaAngDt;
                    sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXn1(:,i)' * this.body1.deltaAngDt + this.contactFrame(:,i)' * this.d;
                    bias = sep * this.biasCoefficient;
                    normalVel = this.contactFrame(:,i) .* this.body1.v + this.body1.w .* this.raXn1(:,i);
                    dlambdaTan(i-1) =  (bias / (this.w1(i)) - sum(normalVel) / (this.w1(i)))*0.8;
                    %dlambdaTan(i-1) = - sum(normalVel) / this.w1(1) * 0.8;
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
                this.body1.v = this.body1.v +  this.delLinVel1 * dlambdas;
                this.body1.w = this.body1.w +  this.angDelta1 * dlambdas;
                this.dlambdaSP = this.dlambdaSP + dlambdas;
            end
        end

        %%
        function applyLambdaSP(this)
                this.body2.v = this.body2.v - this.delLinVel2 * this.dlambdaSP;
                this.body2.w = this.body2.w - this.angDelta2 * this.dlambdaSP;
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
