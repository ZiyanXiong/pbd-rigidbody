classdef JointHinge2 < Joint
    properties
        axis % Axis in the local coordiante of body 1
        xl1  % Joint position in the local coordiante of body1
        torques
        limits
        limitSigns
    end

    methods
        function this = JointHinge2(body1, body2, ground, xl1, axis, torques)
            this = this@Joint(body1,body2, ground);
            this.xl1 = xl1;
            this.axis = axis / norm(axis);
            this.constraintNum = 2;
            this.lambdaLen = 6;
            this.constraints{end+1} = apbd.ConFix(this.body1,this.body2, xl1, this.axis);
            this.constraints{end+1} = apbd.ConRotate(this.body1,this.body2, xl1, this.axis);
            this.limits = [pi,-pi];
            this.limitSigns = [];

             if(nargin < 6)
                 this.torques = [];
             else
                 this.torques = torques;
             end
        end

        function init(this,h,hs,~)
            n = this.lambdaLen;
            this.J1I = zeros(n,6);
            this.J2I = zeros(n,6);
            this.b = zeros(n,1);
            for i = 1:this.constraintNum
                this.constraints{i}.init(h,hs,this.limits(1),this.limits(2),true);
            end
            this.limitSigns = [0;0;0;this.constraints{2}.limitSign;0;0;];
        end

        %%
        function setLimits(this, limitHight, limitLow)
            this.limits = [limitHight, limitLow];
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
                    rows = 3*(i-1) + 1: 3*i;
                    if(i~=1)
                        this.J1I(rows,1:3) = this.constraints{i}.raXnI';
                        this.J1I(rows,4:6) = this.constraints{i}.contactFrame' ./ sqrt(this.body1.Mp);
                    else
                        this.J1I(rows,1:3) = this.constraints{i}.nI1';
                        this.J1I(rows,4:6) = zeros(3,3);
                    end
                    this.b(rows) = -this.constraints{i}.evalCs();
                end
            else
                for i = 1:this.constraintNum
                    if(i==1)
                        rows = 1:3;
                        this.J1I(rows,1:3) = this.constraints{i}.raXnI1';
                        this.J1I(rows,4:6) = this.constraints{i}.contactFrame' ./ sqrt(this.constraints{i}.body1.Mp);
    
                        this.J2I(rows,1:3) = -this.constraints{i}.raXnI2';
                        this.J2I(rows,4:6) = -this.constraints{i}.contactFrame' ./ sqrt(this.constraints{i}.body2.Mp);
                        Cs = this.constraints{i}.evalCs();
                        this.b(rows) = -Cs;
                    else
                        rows = 4:6;
                        this.J1I(rows,1:3) = this.constraints{i}.nI1';
                        this.J1I(rows,4:6) = zeros(3,3);

                        this.J2I(rows,1:3) = -this.constraints{i}.nI2';
                        this.J2I(rows,4:6) = zeros(3,3);
                        Cs = this.constraints{i}.evalCs();
                        this.b(rows) = -Cs;
                    end
                end
            end
        end

        %%
        function compute_b(this)
            for i = 1:this.constraintNum
                if(i==1)
                    rows = 1:3;
                    Cs = this.constraints{i}.evalCs();
                    this.b(rows) = -Cs;
                else
                    rows = 4:5;
                    Cs = this.constraints{i}.evalCs();
                    this.b(rows) = -Cs(2:3);
                end
            end
        end

        %%
        function compute_d(this)
            for i = 1:this.constraintNum
                if(i == 1)
                    rows = 1:3;
                    this.d(rows) = this.constraints{i}.contactFrame'* this.constraints{i}.dt;
                else
                    rows = 4:6;
                    this.d(rows) = this.constraints{i}.dt;
                end
            end
        end

        %%
        function applyLambdas(this, lambdas)
            for i = 1:this.constraintNum
                if(i == 1)
                    rows = 1:3;
                    this.constraints{i}.applyLambda(lambdas(rows));
                else
                    rows = 4:6;
                    this.constraints{i}.applyLambda(lambdas(rows));
                end
            end
        end

        %%
        function recordTorques(this, timestep)
            this.torques(timestep) = this.constraints{2}.lambda(1);
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