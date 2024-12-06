classdef Collision < handle
    properties
        body1   
        body2   % If this is ground collision, body2 will be null
        contacts
        contactNum
        constraints
        ground  %If this is a ground collision
        broken  %If we need to do collision detecitno again
        index   % Global begining index for each collision 
        mIndces % Indices in the matrix
        nextColl % List of next collisions
        layer    % Layer of this collision
        J1I
        J2I     % If this is ground collision, J2 will be null
        b
    end

    methods
        function this = Collision(body1, body2, ground)
            this.contactNum = 0;
            this.contacts = {Contact(), Contact(), Contact(), Contact(), Contact(), Contact(), Contact(), Contact()};
            this.constraints = {};
            this.body1 = body1;
            this.body2 = body2;
            this.ground = ground;
            this.broken = true;
            this.index = 0;
        end

        %%
        function setContacts(this,cdata)
            this.contactNum = length(cdata);
            for i = 1 : this.contactNum
                this.contacts{i}.setData(cdata(i));
            end
            this.J1I = zeros(3*this.contactNum,6);
            this.J2I = zeros(3*this.contactNum,6);
            this.b = zeros(3*this.contactNum,1);
        end

        %%
        function getConstraints(this)
            this.constraints = {};
            for i = 1 : this.contactNum
                if this.ground
                    this.constraints{end+1} = apbd.ConCollGroundRigid(this.body1,this.contacts{i}, this);
                else
                    this.constraints{end+1} = apbd.ConCollRigidRigid(this.body1, this.body2, this.contacts{i}, this);
                end
            end
        end

        %%
        function computeJ_b(this)
            if(this.ground)
                %I1sqrt = 1 ./ sqrt([this.body1.Mr; ones(3,1)*this.body1.Mp]);
                for i = 1:this.contactNum
                    rows = 3*(i-1) + 1: 3*i;
                    this.J1I(rows,1:3) = this.constraints{i}.raXnI';
                    this.J1I(rows,4:6) = this.constraints{i}.contactFrame' ./ sqrt(this.body1.Mp);
                    this.b(rows) = -this.constraints{i}.evalCs();
                end
            else
                for i = 1:this.contactNum
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
        function solveCollisionNor(this, minPenetration, withSP)
            for i = 1 : this.contactNum
                this.constraints{i}.solveNorPos(minPenetration, withSP);
            end
        end

        %%
        function solveCollisionTan(this, withSP)
            for i = 1 : this.contactNum
                this.constraints{i}.solveTanVel(withSP);
            end
        end

        %%
        function applyLambdaSP(this)
            for i = 1 : this.contactNum
                this.constraints{i}.applyLambdaSP();
            end
        end

        %%
        function initConstraints(this, hs, biasCoeff)
            for i = 1 : this.contactNum
                this.constraints{i}.init(hs, biasCoeff);
            end
            if(this.body1.layer < this.body2.layer)
                temp = this.body1;
                this.body1 = this.body2;
                this.body2 = temp;
            end
        end

        %%
        function draw(this)
            for i = 1 : this.contactNum
                this.constraints{i}.draw();
            end
        end
    end
end