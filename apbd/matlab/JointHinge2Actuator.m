classdef JointHinge2Actuator < Joint
    properties
        torques
        targetsTheta
        targetsW
    end

    methods
        function this = JointHinge2Actuator(body1, body2, ground, xl1, axis, targets, targetsW, kp, kd)
            this = this@Joint(body1,body2, ground);
            this.constraintNum = 2;
            this.lambdaLen = 6;
            axis = axis / norm(axis);
            this.constraints{end+1} = apbd.ConFix(this.body1,this.body2, xl1, axis);
            this.constraints{end+1} = apbd.ConRotate(this.body1,this.body2, xl1, axis, kp, kd);
            this.targetsTheta = targets;
            this.targetsW = targetsW;
            this.torques = zeros(length(this.targetsTheta),1);
        end

        function init(this,h,hs,timestep)
            n = this.lambdaLen;
            this.J1I = zeros(n,6);
            this.J2I = zeros(n,6);
            this.b = zeros(n,1);
            for i = 1:this.constraintNum
                this.constraints{i}.init(h,hs,this.targetsTheta(timestep), this.targetsW(timestep));
            end
        end

        %%
        function [f,t]= applyForceTorque(this,f,t,timestep)
            if(~isempty(this.torques))
                f(this.body1.index,:) = f(this.body1.index,:) + zeros(1,3);
                f(this.body2.index,:) = f(this.body2.index,:) + zeros(1,3);
                t(this.body1.index,:) = t(this.body1.index,:) + this.torques(timestep,:);
                t(this.body2.index,:) = t(this.body2.index,:) - this.torques(timestep,:);
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
                        this.J1I(rows,4:6) = zeros;

                        this.J2I(rows,1:3) = -this.constraints{i}.nI2';
                        this.J2I(rows,4:6) = zeros;
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
                    rows = 4:6;
                    Cs = this.constraints{i}.evalCs();
                    this.b(rows) = -Cs;
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