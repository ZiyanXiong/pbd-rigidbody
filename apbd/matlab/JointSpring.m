classdef JointSpring < Joint
    properties
        stiffness
        torques
    end

    methods
        function this = JointSpring(body1, body2, ground, xl1, xl2, dis, stiffness)
            this = this@Joint(body1,body2, ground);
            this.stiffness = stiffness;
            this.constraintNum = 1;
            this.lambdaLen = 1;
            this.constraints{end+1} = apbd.ConFixDistance(body1,body2, xl1, xl2, dis);
            this.torques = [];
        end

        function init(this,h,hs,~)
            n = 1;
            this.J1I = zeros(n,6);
            this.J2I = zeros(n,6);
            this.b = zeros(n,1);
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
            if(this.ground)
                %I1sqrt = 1 ./ sqrt([this.body1.Mr; ones(3,1)*this.body1.Mp]);
                for i = 1:this.constraintNum
                    rows = 1;
                    this.J1I(rows,1:3) = this.constraints{i}.raXnI';
                    this.J1I(rows,4:6) = this.constraints{i}.contactFrame' ./ sqrt(this.body1.Mp);
                    this.b(rows) = -this.constraints{i}.evalCs();
                end
            else
                for i = 1:this.constraintNum
                    rows = 1;
                    this.J1I(rows,1:3) = this.constraints{i}.raXnI1';
                    this.J1I(rows,4:6) = this.constraints{i}.contactFrame' ./ sqrt(this.constraints{i}.body1.Mp);

                    this.J2I(rows,1:3) = -this.constraints{i}.raXnI2';
                    this.J2I(rows,4:6) = -this.constraints{i}.contactFrame' ./ sqrt(this.constraints{i}.body2.Mp);
                    this.b(rows) = -this.constraints{i}.evalCs();
                end
            end
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