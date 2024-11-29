classdef Collision < handle
    properties
        body1   
        body2   % If this is ground collision, body2 will be null
        contacts
        contactNum
        constraints
        ground  %If this is a ground collision
        broken  %If we need to do collision detecitno again
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
        end

        %%
        function setContacts(this,cdata)
            this.contactNum = length(cdata);
            for i = 1 : this.contactNum
                this.contacts{i}.setData(cdata(i));
            end
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
        end

        %%
        function draw(this)
            for i = 1 : this.contactNum
                this.constraints{i}.draw();
            end
        end
    end
end