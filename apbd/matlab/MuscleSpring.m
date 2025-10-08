classdef MuscleSpring < Muscle
    properties
        stiffness
        numBodies
        torques
    end

    methods
        function this = MuscleSpring(bodies, ground, points, dis, stiffness)
            this = this@Muscle(bodies, ground);
            this.stiffness = stiffness;
            this.constraintNum = 1;
            this.lambdaLen = 1;
            numPoints = length(points);
            this.numBodies = length(bodies);
            viaPoints = points(:,2:end-1);
            dis = dis - norm(points(:,1)-points(:,2)) - norm(points(:,numPoints - 1)- points(:,numPoints));
            this.constraints{end+1} = apbd.ConFixDistanceViaPoint(bodies, viaPoints, dis);
            this.torques = [];
        end

        function init(this,h,hs,~)
            this.JIs = zeros(this.viaPointNum, 6);
            this.b = zeros(1,1);
            this.compliance = 1/ this.stiffness / (h*h);
            for i = 1:this.constraintNum
                this.constraints{i}.init(h,hs);
            end
        end

        %%
        function [f,t]= applyForceTorque(this,f,t,timestep)
            if(~isempty(this.torques))
                f(this.body1.index,:) = f(this.body1.index,:) + zeros(1,3);
                f(this.body2.index,:) = f(this.body2.index,:) + zeros(1,3);
                axisW = this.constraints{2}.body1.transformVector(this.axis);
                t(this.body1.index,:) = t(this.body1.index,:) + this.torques(timestep) * axisW';
                t(this.body2.index,:) = t(this.body2.index,:) - this.torques(timestep) * axisW';
            end
        end

        %%
        function computeJ_b(this)
            for i = 1:this.constraintNum
                for j = 1:this.numBodies
                    this.JIs(j,1:3) = this.constraints{i}.raXnIs(:,j)';
                    this.JIs(j,4:6) = this.constraints{i}.nw(:,j)' ./ sqrt(this.constraints{i}.bodies{j}.Mp);
                end
            end
            this.b(1) =  -this.constraints{1}.evalCs();
        end

        %%
        function compute_b(this)
            for i = 1:this.constraintNum
                rows = 1;
                this.b(rows) = -this.constraints{i}.evalCs();
            end
        end

        %%
        function compute_d(this)
            for i = 1:this.constraintNum
                rows = 1;
                %this.d(rows) = this.constraints{i}.contactFrame'* this.constraints{i}.dt;
                this.d(rows) = 0;
            end
        end

        %%
        function applyLambdas(this, lambdas)
            for i = 1:this.constraintNum
                rows = 1;
                this.constraints{i}.applyLambda(lambdas(rows));
            end
        end

        %%
        function recordTorques(this, ~)
            this.torques = [];
        end

        %%
        function solveCollisionNor(this, withSP)
            for i = 1 : this.constraintNum
                this.constraints{i}.solveNorPos(withSP);
            end
        end

        %%
        function solveCollisionTan(this,withSP)
            for i = 1 : this.constraintNum
                this.constraints{i}.solveTanPos(withSP);
            end
        end
    end
end