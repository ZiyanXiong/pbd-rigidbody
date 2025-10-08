classdef ConFixDistanceViaPoint < apbd.ConColl
	%ConCollRigidRigid Collision between two rigid bodies

	properties
        bodies
        viaPoints
        numBodies
        dis % Distance between 2 points
        dt % d/h

        contactFrame
        ws   % Generalized mass vector (3x1)
        delLinVels % Unit change for linear velocity matrix(3x3)
        angDeltas % Unit change for angular velocity matrix(3x3)
        raXns  % ra X nw * sqrt(I^(-1)) matrix(3x3)
        raXnIs  % ra X nw * sqrt(I^(-1)) matrix(3x3)
        dlambdas

        biasCoefficient
	end

	methods
		%%
        function this = ConFixDistanceViaPoint(bodies, viaPoints, dis)
			this = this@apbd.ConColl();
			this.bodies = bodies;
            this.numBodies = length(bodies);
            this.viaPoints = viaPoints;
			this.nw = zeros(3, this.numBodies);
            this.dis = dis;

            this.ws = zeros(1, this.numBodies);
            this.raXns = zeros(3, this.numBodies);
            this.raXnIs = zeros(3, this.numBodies);
            this.delLinVels = zeros(3, this.numBodies);
            this.angDeltas = zeros(3, this.numBodies);

             this.dlambdas = zeros(1,1);
		end

		%%
		function init(this,h,hs) 
            this.d = - this.dis;
            this.nw = zeros(3, this.numBodies);
            for i = 1:this.numBodies
                if(i < this.numBodies)
                    d = this.bodies{i}.transformPoint(this.viaPoints(:,i)) - this.bodies{i+1}.transformPoint(this.viaPoints(:,i+1));
                    dNorm = norm(d);
                    this.nw(:,i) = this.nw(:,i) + d / dNorm;
                    this.nw(:,i+1) = this.nw(:,i+1) - d / dNorm;
                    this.d = this.d + dNorm;
                end
    			m = this.bodies{i}.Mp;
			    I = this.bodies{i}.Mr;
			    q = this.bodies{i}.x(1:4);
			    rl = this.viaPoints(:,i);
                nl = se3.qRotInv(q, this.nw(:,i));
			    rnl = se3.cross(rl,nl);
                this.raXnIs(:,i) = se3.qRot(q,(sqrt(I).\rnl));
                this.raXns(:,i) = se3.qRot(q,rnl);
			    this.ws(i) = (1/m)*(nl'*nl) + this.raXnIs(:,i)' * this.raXnIs(:,i);
                this.delLinVels(:,i) = this.nw(:,i) / m;
                this.angDeltas(:,i) = se3.qRot(q,(I.\rnl));
            end
            this.dt = this.d / h;
            this.biasCoefficient = -1 / hs;
            this.lambda = zeros(1,1);
        end

        %%
        function layer = getLayer(this)
            layer = this.body1.layer + this.body2.layer;
        end

        %%
        function Cs = evalCs(this)
            Cs = this.dt;
            for i = 1:this.numBodies
                Cs = Cs + this.nw(:,i)' * this.bodies{i}.v + this.raXns(:,i)' * this.bodies{i}.w;
            end
        end

        %%
        function applyLambda(this, dlambdas)
            this.lambda = this.lambda + dlambdas;
            for i = 1:this.numBodies
                this.bodies{i}.v = this.bodies{i}.v +  this.delLinVels(:,i) * dlambdas;
                this.bodies{i}.w = this.bodies{i}.w +  this.angDeltas(:,i) * dlambdas;
            end 
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
            xs = zeros(3, this.numBodies);
            for i = 1:this.numBodies
			    xs(:,i) = this.bodies{i}.transformPoint(this.viaPoints(:,i));
            end
            for i = 1:this.numBodies
			    plot3(xs(1,i),xs(2,i),xs(3,i),'go');
                if(i< this.numBodies)
                    plot3(xs(1,i:i+1),xs(2,i:i+1),xs(3,i:i+1),'g');
                end
            end
		end
	end
end
