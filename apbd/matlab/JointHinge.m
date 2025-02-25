classdef JointHinge < Joint
    properties
        axis
        xw
    end

    methods
        function this = JointHinge(body1, body2, ground,xw,axis)
            this = this@Joint(body1,body2, ground);
            this.xw = xw;
            this.axis = axis;
            this.constraintNum = 1;
            xw1 = this.xw + this.axis*0;
            %xw2 = this.xw - this.axis*4;
            this.constraints{end+1} = apbd.ConFix(this.body1,this.body2, xw1, [0 0 1]');
            %this.constraints{end+1} = apbd.ConFix(this.body1,this.body2, xw2, this.axis);
        end

        function init(this,h,hs)
            n = 3*this.constraintNum;
            this.J1I = zeros(n,6);
            this.J2I = zeros(n,6);
            this.b = zeros(n,1);
            for i = 1:this.constraintNum
                this.constraints{i}.init(h,hs);
            end
        end

        %%
        function computeJ_b(this)
            if(this.ground)
                %I1sqrt = 1 ./ sqrt([this.body1.Mr; ones(3,1)*this.body1.Mp]);
                for i = 1:this.constraintNum
                    rows = 3*(i-1) + 1: 3*i;
                    this.J1I(rows,1:3) = this.constraints{i}.raXnI';
                    this.J1I(rows,4:6) = this.constraints{i}.contactFrame' ./ sqrt(this.body1.Mp);
                    this.b(rows) = -this.constraints{i}.evalCs();
                end
            else
                for i = 1:this.constraintNum
                    rows = 3*(i-1) + 1: 3*i;
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
                rows = 3*(i-1) + 1: 3*i;
                this.b(rows) = -this.constraints{i}.evalCs();
            end
        end

        %%
        function compute_d(this)
            for i = 1:this.constraintNum
                rows = 3*(i-1) + 1: 3*i;
                this.d(rows) = this.constraints{i}.contactFrame'* this.constraints{i}.dt;
            end
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