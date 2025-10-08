classdef ConFix < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
		body1
		body2
		x1 % Position wrt body 1 (3x1)
		x2 % Position wrt body 2 (3x1)
        nl % Normal wrt body1 (3x1)
        dt % d/h

        contactFrame
        w1   % Generalized mass vector (3x1)
        delLinVel1 % Unit change for linear velocity matrix(3x3)
        angDelta1 % Unit change for angular velocity matrix(3x3)
        raXn1  % ra X nw * sqrt(I^(-1)) matrix(3x3)
        raXnI1  % ra X nw * sqrt(I^(-1)) matrix(3x3)

        w2   % Generalized mass vector (3x1)
        delLinVel2 % Unit change for linear velocity matrix(3x3)
        angDelta2 % Unit change for angular velocity matrix(3x3)
        raXn2  % ra X nw * sqrt(I^(-1)) matrix(3x3)
        raXnI2  % ra X nw * sqrt(I^(-1)) matrix(3x3)
        dlambdas

        biasCoefficient
	end

	methods
		%%
        function this = ConFix(body1,body2, xl1, nl)
			this = this@apbd.ConColl();
			this.body1 = body1;
			this.body2 = body2;
			this.nl = nl;
			this.x1 = xl1;
			this.x2 = body2.invTransformPoint(body1.transformPoint(this.x1));

            this.contactFrame = zeros(3,3);

            this.w1 = zeros(3,1);
            this.raXn1 = zeros(3,3);
            this.delLinVel1 = zeros(3,3);
            this.angDelta1 = zeros(3,3);

            this.w2 = zeros(3,1);
            this.raXn2 = zeros(3,3);
            this.delLinVel2 = zeros(3,3);
            this.angDelta2 = zeros(3,3);

             this.dlambdas = zeros(3,1);
		end

		%%
		function init(this,h,hs,~,~,~) 
            this.d = this.body1.transformPoint(this.x1) - this.body2.transformPoint(this.x2);
            this.dt = this.d / h;
            this.biasCoefficient = -1 / hs;

            this.lambda = zeros(3,1);
            this.nw = this.body1.transformVector(this.nl);
            [tanx,tany] = apbd.ConColl.generateTangents(this.nw);
            this.contactFrame = [this.nw, tanx, tany];

			m1 = this.body1.Mp;
			I1 = this.body1.Mr;
			q1 = this.body1.x(1:4);
			rl1 = this.x1;

			m2 = this.body2.Mp;
			I2 = this.body2.Mr;
			q2 = this.body2.x(1:4);
			rl2 = this.x2;
            
            for i = 1:3
                nl1 = se3.qRotInv(q1, this.contactFrame(:,i));
			    rnl1 = se3.cross(rl1,nl1);
                this.raXnI1(:,i) = se3.qRot(q1,(sqrt(I1).\rnl1));
                this.raXn1(:,i) = se3.qRot(q1,rnl1);
			    this.w1(i) = (1/m1) + this.raXnI1(:,i)' * this.raXnI1(:,i);
                this.delLinVel1(:,i) = this.contactFrame(:,i) / m1;
                this.angDelta1(:,i) = se3.qRot(q1,(I1.\rnl1));

                nl2 = se3.qRotInv(q2, this.contactFrame(:,i));
                rnl2 = se3.cross(rl2,nl2);
                this.raXnI2(:,i) = se3.qRot(q2,(sqrt(I2).\rnl2));
                this.raXn2(:,i) = se3.qRot(q2,rnl2);
			    this.w2(i) = (1/m2) + this.raXnI2(:,i)' * this.raXnI2(:,i);
                this.delLinVel2(:,i) = this.contactFrame(:,i) / m2;
                this.angDelta2(:,i) = se3.qRot(q2,(I2.\rnl2));
            end
        end

        %%
        function layer = getLayer(this)
            layer = this.body1.layer + this.body2.layer;
        end

        %%
        function Cs = evalCs(this)
            Cs = this.contactFrame' * (this.body1.v - this.body2.v) + this.raXn1' * this.body1.w - this.raXn2' * this.body2.w + this.contactFrame'* this.dt;
        end

        %%
        function applyLambda(this, dlambdas)
            this.lambda = this.lambda + dlambdas;
            this.body1.v = this.body1.v +  this.delLinVel1 * dlambdas;
            this.body1.w = this.body1.w +  this.angDelta1 * dlambdas;
            this.body2.v = this.body2.v - this.delLinVel2 * dlambdas;
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
			x = this.body1.transformPoint(this.x1);
			plot3(x(1),x(2),x(3),'ro');
			x = this.body2.transformPoint(this.x2);
			plot3(x(1),x(2),x(3),'go');
		end
	end
end
