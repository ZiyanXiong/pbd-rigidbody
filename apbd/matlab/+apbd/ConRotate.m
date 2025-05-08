classdef ConRotate < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
		body1
		body2
        x1 % Position wrt body 1 (3x1)
		cf1 % Contact frame wrt body 1 (3x1)
		cf2 % Contact frame wrt body 2 (3x1)
        dt % d/h
        kp
        kd

        contactFrame
        w1   % Generalized mass vector (3x1)
        angDelta1 % Unit change for angular velocity matrix(3x3)
        nI1  % nw * sqrt(I^(-1)) matrix(3x3)

        w2   % Generalized mass vector (3x1)
        angDelta2 % Unit change for angular velocity matrix(3x3)
        nI2  % nw * sqrt(I^(-1)) matrix(3x3)
        dlambdas

        biasCoefficient
	end

	methods
		%%
        function this = ConRotate(body1, body2, x, axis, kp, kd)
            if nargin < 5
                kp = 1;
                kd =0;
            end

			this = this@apbd.ConColl();
			this.body1 = body1;
			this.body2 = body2;
            this.x1 = x;
            [tanx,tany] = apbd.ConColl.generateTangents(axis);
            this.contactFrame = [axis, tanx, tany];
			this.cf1 = this.contactFrame;
            for i = 1:3
    			this.cf2(:,i) = body2.invTransformVector(body1.transformVector(this.cf1(:,i)));
            end
            
            this.w1 = zeros(1,1);
            this.nI1 = zeros(3,3);
            this.angDelta1 = zeros(3,3);

            this.w2 = zeros(1,1);
            this.nI2 = zeros(3,3);
            this.angDelta2 = zeros(3,3);

            this.dlambdas = zeros(3,1);

            this.kp = kp;
            this.kd = kd;
		end

		%%
		function init(this,h,hs,thetaTarget, wTarget)
			I1 = this.body1.Mr;
			q1 = this.body1.x(1:4);

			I2 = this.body2.Mr;
			q2 = this.body2.x(1:4);

            for i = 1:3
                nl1 = this.cf1(:,i);
                this.nI1(:,i) = se3.qRot(q1, sqrt(I1).\nl1);
                this.angDelta1(:,i) = se3.qRot(q1,(I1.\ nl1));
                
                nl2 = se3.qRotInv(q2, this.contactFrame(:,i));
                this.nI2(:,i) = se3.qRot(q2,sqrt(I2).\nl2);
                this.angDelta2(:,i) = se3.qRot(q2,(I2.\ nl2));

                this.contactFrame(:,i) = this.body1.transformVector(this.cf1(:,i));
            end
          
            dqAlign = se3.computeDq(this.body2.transformVector(this.cf2(:,1)), this.body1.transformVector(this.cf1(:,1)));
            dqTarget = se3.computeDq(se3.qRot(dqAlign,this.body2.transformVector(this.cf2(:,2))), this.body1.transformVector(this.cf1(:,2)));

            dtheta = se3.dqToDeltaTheta(dqAlign);
            this.dt = this.contactFrame' * dtheta / h;
            dthetaTarget = this.contactFrame(:,1)'*se3.dqToDeltaTheta(dqTarget);
            
            angVelocity = this.contactFrame(:,1)'*(this.body1.w - this.body2.w);
            a = (h/(h*(h*this.kp+this.kd)*(this.contactFrame(:,1)'* (this.angDelta1(:,1) - this.angDelta2(:,1))) + 1));
            %a = (h/(h*(h*this.kp+this.kd) + 1));
            this.dt(1) = a * (this.kp*((dthetaTarget - thetaTarget) + angVelocity * h) + ...
                this.kd*(angVelocity - wTarget)) - angVelocity;            
            this.biasCoefficient = -1 / hs;

            this.lambda = zeros(3,1);
        end

        %%
        function layer = getLayer(this)
            layer = this.body1.layer + this.body2.layer;
        end

        %%
        function Cs = evalCs(this)
            Cs =  this.contactFrame' * (this.body1.w - this.body2.w) + this.dt;
        end

        %%
        function applyLambda(this, dlambdas)
            this.lambda = this.lambda + dlambdas;
            this.body1.w = this.body1.w +  this.angDelta1 * dlambdas;
            this.body2.w = this.body2.w - this.angDelta2 * dlambdas;
        end

		%%
        function solveNorPos(this, withSP)
            if(~withSP)
                sep = this.nw' * this.body1.deltaLinDt + this.raXn1(:,1)' * this.body1.deltaAngDt + this.nw' * this.d;
                sep = sep - (this.nw' * this.body2.deltaLinDt + this.raXn2(:,1)' * this.body2.deltaAngDt);
                bias = sep * this.biasCoefficient;
                %normalVel = this.nw' * this.body.computePointVel(this.xl);
                normalVel = this.nw .* this.body1.v + this.body1.w .* this.raXn1(:,1);
                normalVel = normalVel - (this.nw .* this.body2.v + this.body2.w .* this.raXn2(:,1));
                this.dlambdas(1) =  bias / (this.w1(1) + this.w2(1)) - sum(normalVel) / (this.w1(1) + this.w2(1));
                this.lambda(1) = this.lambda(1) + this.dlambdas(1);
                this.body1.v = this.body1.v + this.dlambdas(1) * this.delLinVel1(:,1);
                this.body1.w = this.body1.w + this.dlambdas(1) * this.angDelta1(:,1);
                this.body2.v = this.body2.v - this.dlambdas(1) * this.delLinVel2(:,1);
                this.body2.w = this.body2.w - this.dlambdas(1) * this.angDelta2(:,1);
            else
                sep = this.nw' * this.body1.deltaLinDt + this.raXn1(:,1)' * this.body1.deltaAngDt + this.nw' * this.dt;
                sep = sep - (this.nw' * this.body2.deltaLinDt + this.raXn2(:,1)' * this.body2.deltaAngDt);
                bias = -sep ;
                %normalVel = this.nw' * this.body.computePointVel(this.xl);
                normalVel = this.nw .* this.body1.v + this.body1.w .* this.raXn1(:,1);
                normalVel = normalVel - (this.nw .* this.body2.v + this.body2.w .* this.raXn2(:,1));
                if(this.body1.layer == this.body2.layer)
                    this.dlambdas(1) =  bias / (this.w1(1)+this.w2(1)) - sum(normalVel) / (this.w1(1)+this.w2(1));
                else
                    this.dlambdas(1) =  bias / this.w1(1) - sum(normalVel) / this.w1(1);
                end
                lambda = this.lambda(1) + this.dlambdas(1);
                if(lambda < 0)
                    this.dlambdas(1) = - this.lambda(1);
                    this.collision.broken = true;
                end
                this.lambda(1) = this.lambda(1) + this.dlambdas(1);
                if(this.body1.layer == this.body2.layer)
                    this.body1.v = this.body1.v + this.dlambdas(1) * this.delLinVel1(:,1);
                    this.body1.w = this.body1.w + this.dlambdas(1) * this.angDelta1(:,1);
                    this.body2.v = this.body2.v - this.dlambdas(1) * this.delLinVel2(:,1);
                    this.body2.w = this.body2.w - this.dlambdas(1) * this.angDelta2(:,1);
                else
                    this.body1.v = this.body1.v + this.dlambdas(1) * this.delLinVel1(:,1);
                    this.body1.w = this.body1.w + this.dlambdas(1) * this.angDelta1(:,1);
                end
            end
        end

		%%
        function solveTanPos(this, withSP)
            if(~withSP)
                dlambdaTan = zeros(2,1);
                for i = 2:3
                    %sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXnI1(:,i)' * this.body1.deltaAngDt;
                    sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXn1(:,i)' * this.body1.deltaAngDt + this.contactFrame(:,i)' * this.d;
                    sep = sep - (this.contactFrame(:,i)' * this.body2.deltaLinDt + this.raXn2(:,i)' * this.body2.deltaAngDt);
                    bias = sep * this.biasCoefficient;
                    normalVel = this.contactFrame(:,i) .* this.body1.v + this.body1.w .* this.raXn1(:,i);
                    normalVel = normalVel - (this.contactFrame(:,i) .* this.body2.v + this.body2.w .* this.raXn2(:,i));
                    dlambdaTan(i-1) =  (bias / (this.w1(i) + this.w2(i)) - sum(normalVel) / (this.w1(i) + this.w2(i)));
                end
                dlambdaTan = [0;dlambdaTan];
                this.lambda = this.lambda + dlambdaTan;
                this.body1.v = this.body1.v +  this.delLinVel1 * dlambdaTan;
                this.body1.w = this.body1.w +  this.angDelta1 * dlambdaTan;
                this.body2.v = this.body2.v - this.delLinVel2 * dlambdaTan;
                this.body2.w = this.body2.w - this.angDelta2 * dlambdaTan;
                this.dlambdas(2:3) = dlambdaTan(2:3);
            else
                dlambdaTan = zeros(2,1);
                for i = 2:3
                    %sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXnI1(:,i)' * this.body1.deltaAngDt;
                    sep = this.contactFrame(:,i)' * this.body1.deltaLinDt + this.raXn1(:,i)' * this.body1.deltaAngDt + this.contactFrame(:,i)' * this.dt;
                    sep = sep - (this.contactFrame(:,i)' * this.body2.deltaLinDt + this.raXn2(:,i)' * this.body2.deltaAngDt);
                    bias = -sep;
                    normalVel = this.contactFrame(:,i) .* this.body1.v + this.body1.w .* this.raXn1(:,i);
                    normalVel = normalVel - (this.contactFrame(:,i) .* this.body2.v + this.body2.w .* this.raXn2(:,i));
                    if(this.body1.layer == this.body2.layer)
                        dlambdaTan(i-1) =  (bias / (this.w1(i) + this.w2(i)) - sum(normalVel) / (this.w1(i) + this.w2(i)));
                    else
                        dlambdaTan(i-1) =  (bias / this.w1(i) - sum(normalVel) / this.w1(i));
                    end
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
                if(this.body1.layer == this.body2.layer)
                    this.body1.v = this.body1.v +  this.delLinVel1 * dlambdaTan;
                    this.body1.w = this.body1.w +  this.angDelta1 * dlambdaTan;
                    this.body2.v = this.body2.v - this.delLinVel2 * dlambdaTan;
                    this.body2.w = this.body2.w - this.angDelta2 * dlambdaTan;
                else
                    this.body1.v = this.body1.v +  this.delLinVel1 * dlambdaTan;
                    this.body1.w = this.body1.w +  this.angDelta1 * dlambdaTan;
                end
                this.dlambdas(2:3) = dlambdaTan(2:3);
            end
        end

        %%
        function applyLambdaSP(this)
            if(this.body1.layer ~= this.body2.layer)
                this.body2.v = this.body2.v - this.delLinVel2 * this.lambda;
                this.body2.w = this.body2.w - this.angDelta2 * this.lambda;
            end
        end

		%%
		function draw(this)
            E(1:3,4) = this.body1.transformPoint(this.x1);
            E(1:3,1:3) = this.contactFrame;
            se3.drawAxis(E);
		end
	end
end
