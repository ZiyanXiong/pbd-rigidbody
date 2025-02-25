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
        b0
        bPrev
        d
        mu
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
            if(this.ground)
                this.mu = this.body1.mu;
            else
                this.mu = 0.5 * (this.body1.mu + this.body2.mu);
            end
        end

        %%
        function setContacts(this,cdata)
            this.contactNum = min(length(cdata),8);
            for i = 1 : this.contactNum
                this.contacts{i}.setData(cdata(i));
            end
            n = 3*this.contactNum;
            this.J1I = zeros(n,6);
            this.J2I = zeros(n,6);
            this.b = zeros(n,1);
            this.b0 = zeros(n,1);
            this.bPrev = zeros(n,1);
            this.d = zeros(n,1);
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
        function compute_b(this)
            for i = 1:this.contactNum
                rows = 3*(i-1) + 1: 3*i;
                this.b(rows) = -this.constraints{i}.evalCs();
            end
        end

        %%
        function compute_b0(this)
            for i = 1:this.contactNum
                rows = 3*(i-1) + 1: 3*i;
                this.b0(rows) = -this.constraints{i}.evalCs();
            end
        end
        %%
        function compute_d(this)
            for i = 1:this.contactNum
                rows = 3*(i-1) + 1: 3*i;
                this.d(rows) = this.constraints{i}.contactFrame'* this.constraints{i}.dt;
            end
        end

        %%
        function solveCollisionNor(this, withSP)
            for i = 1 : this.contactNum
                this.constraints{i}.solveNorPos(withSP);
            end
        end

        %%
        function solveCollisionTan(this,withSP)
            for i = 1 : this.contactNum
                this.constraints{i}.solveTanPos(withSP);
            end
        end

        %%
        function solveCollisionNorVel(this, substeps)
            for i = 1 : this.contactNum
                this.constraints{i}.solveNorVel(substeps);
            end
        end

        %%
        function solveCollisionTanVel(this, substeps)
            for i = 1 : this.contactNum
                this.constraints{i}.solveTanVel(substeps);
            end
        end

        %%
        function iscon = isConverged(this)
            this.bPrev = this.b;
            this.compute_b();
            iscon = norm(this.b-this.bPrev) < 1e-6;
        end

        %%
        function issta = isStable(this)
            deltab = this.b-this.b0;
            bn = this.b(1:3:end);
            issta = all(bn(deltab(1:3:end)<-1e-1) > -1e-1);
        end

        %%
        function resetLambdas(this)
            for i = 1 : this.contactNum
               this.constraints{i}.lambda=zeros(3,1);
               this.constraints{i}.dlambdas=zeros(3,1);
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
        function applyLambdaSP(this)
            for i = 1 : this.contactNum
                this.constraints{i}.applyLambdaSP();
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