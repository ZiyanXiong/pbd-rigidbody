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
        d
        mu

        %% Varibles for GPQP
        lambda
        lambdac
        lambdad
        t_bar
        g
        p
        Ax       % Buffer for storing the resutls of A*x (Matrix dot vector)
        freeIndex  %(111:Free, 100:Degenerate, 000: Fixed)

        %% Varibles for CG
        r_cg
        b_cg
        J1I_cg
        J2I_cg
        Minv_cg
        g_cg
        d_cg
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
            this.d = zeros(n,1);

            this.lambda = zeros(n,1);
            this.lambdac = zeros(n,1);
            this.lambdad = zeros(n,1);

            this.t_bar = zeros(n,1);
            this.g = zeros(n,1);
            this.p = zeros(n,1);
            this.Ax = zeros(n,1);
            this.freeIndex = false(n,1);

            this.J1I = zeros(n,6);
            this.J2I = zeros(n,6);
            this.r_cg = zeros(n,1);
            this.b_cg = zeros(n,1);
            this.Minv_cg = zeros(n,1);
            this.g_cg = zeros(n,1);
            this.d_cg = zeros(n,1);
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
        function compute_d(this)
            for i = 1:this.contactNum
                rows = 3*(i-1) + 1: 3*i;
                this.d(rows) = this.constraints{i}.contactFrame'* this.constraints{i}.dt;
            end
        end

        %%
        function solveCollisionNor(this, minPenetration)
            for i = 1 : this.contactNum
                this.constraints{i}.solveNorPos(minPenetration);
            end
        end

        %%
        function solveCollisionTan(this)
            for i = 1 : this.contactNum
                this.constraints{i}.solveTanVel();
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

        %% Funcitons for GPQP
        %%
        function compute_LTlambda(this)
            this.body1.LTx = this.body1.LTx + this.J1I' * this.lambda;
            this.body2.LTx = this.body2.LTx + this.J2I' * this.lambda;
        end

        %%
        function compute_degenerate_J1I_J2I_b(this)
            this.J1I_cg = this.J1I;
            this.J2I_cg = this.J2I;
            this.b_cg = this.b;
            for i = 1:3:3*this.contactNum
                if(this.freeIndex(i) && ~this.freeIndex(i+1))
                    this.J1I_cg(i,:) = this.lambdad(i:i+2)' * this.J1I(i:i+2,:);
                    this.J2I_cg(i,:) = this.lambdad(i:i+2)' * this.J2I(i:i+2,:);
                    this.b_cg(i) = this.lambdad(i:i+2)' * this.b(i:i+2);
                    this.lambda(i) = this.lambdad(i:i+2)' * this.lambda(i:i+2);
                end
            end
            %this.Minv_cg = ones(3*this.contactNum,1);
            this.Minv_cg = 1 ./ diag(this.J1I_cg * this.J1I_cg' + this.J2I_cg * this.J2I_cg');
        end

        %%
        function compute_degenerate_LTlambda(this)
            this.body1.LTx = this.body1.LTx + this.J1I_cg(this.freeIndex,:)' * this.lambda(this.freeIndex);
            this.body2.LTx = this.body2.LTx + this.J2I_cg(this.freeIndex,:)' * this.lambda(this.freeIndex);
        end

        %%
        function compute_LTd_cg(this)
            this.body1.LTx = this.body1.LTx + this.J1I_cg(this.freeIndex,:)' * this.d_cg(this.freeIndex);
            this.body2.LTx = this.body2.LTx + this.J2I_cg(this.freeIndex,:)' * this.d_cg(this.freeIndex);
        end

        %%
        function compute_LTp(this)
            this.body1.LTx = this.body1.LTx + this.J1I' * this.p;
            this.body2.LTx = this.body2.LTx + this.J2I' * this.p;
        end

        %%
        function compute_LLTx(this)
            this.Ax = this.J1I * this.body1.LTx + this.J2I * this.body2.LTx;
        end

        %%
        function compute_degenerate_LLTx(this)
            this.Ax = this.J1I_cg * this.body1.LTx + this.J2I_cg * this.body2.LTx;
        end

        %%
        function tbar_i = compute_tbar(this)
            for i = 1:3:3*this.contactNum
                t = Collision.rayConeIntersection(this.lambda(i:i+2),-this.g(i:i+2),this.mu);
                this.t_bar(i) = t;
                if((abs(this.g(i)) < 1e-9 && abs(this.lambda(i))<1e-9))
                    this.t_bar(i+1) = 0;
                elseif(this.g(i)>0)
                    this.t_bar(i+1) = this.lambda(i) / this.g(i);
                else
                    this.t_bar(i+1) = Inf;
                end
                this.t_bar(i+2) = Inf;
            end

            tbar_i = this.t_bar;
        end

        %%
        function compute_lambdac(this, t)
            for i = 1:3:3*this.contactNum
                if(t < this.t_bar(i))
                    this.lambdac(i:i+2) = this.lambda(i:i+2) - t*this.g(i:i+2);
                elseif(t < this.t_bar(i+1))
                    gp = this.g(i:i+2);
                    if(norm(this.lambda(i:i+2))<1e-9)
                        this.lambdad(i:i+2) = [1 0 0]';
                    else
                        this.lambdad(i:i+2) = this.lambda(i:i+2) / norm(this.lambda(i:i+2));
                    end
                    gp = (this.lambdad(i:i+2)'*gp)*this.lambdad(i:i+2);
                    this.lambdac(i:i+2) = this.lambda(i:i+2) - this.t_bar(i)*this.g(i:i+2) - (t - this.t_bar(i))*gp;
                else
                    gp = this.g(i:i+2);
                    if(norm(this.lambda(i:i+2))<1e-9)
                        this.lambdad(i:i+2) = [1 0 0]';
                    else
                        this.lambdad(i:i+2) = this.lambda(i:i+2) / norm(this.lambda(i:i+2));
                    end
                    gp = (this.lambdad(i:i+2)'*gp)*this.lambdad(i:i+2);
                    this.lambdac(i:i+2) = this.lambda(i:i+2) -this.t_bar(i)*this.g(i:i+2) + (this.t_bar(i+1) - this.t_bar(i))*gp;
                end
            end
        end

        %%
        function compute_lambdad(this, t)
            for i = 1:3:3*this.contactNum
                if(t<this.t_bar(i))
                    this.freeIndex(i:i+2) = true;
                elseif(t<this.t_bar(i+1))
                    this.freeIndex(i) = true;
                    this.freeIndex(i+1:i+2) = false;
                else
                    this.freeIndex(i:i+2) = false;
                end

                if(norm(this.lambdac(i:i+2))<1e-9)
                    this.lambdad(i:i+2) = [1 0 0]';
                else
                    this.lambdad(i:i+2) = this.lambdac(i:i+2) / norm(this.lambdac(i:i+2));
                end
            end
        end

        %%
        function compute_p(this,t)
            for i = 1:3:3*this.contactNum
                if(t < this.t_bar(i))
                    this.p(i:i+2) = -this.g(i:i+2);
                elseif(t < this.t_bar(i+1))
                    gp = this.g(i:i+2);
                    if(norm(this.lambda(i:i+2))<1e-9)
                        this.lambdad(i:i+2) = [1 0 0]';
                    else
                        this.lambdad(i:i+2) = this.lambda(i:i+2) / norm(this.lambda(i:i+2));
                    end
                    gp = (this.lambdad(i:i+2)'*gp)*this.lambdad(i:i+2);
                    this.p(i:i+2) = -gp;
                else
                    this.p(i:i+2) = zeros(3,1);
                end
            end
        end

        %%
        function project(this)

            for i = 1:3:3*this.contactNum
                if(this.lambda(i) < 0)
                    this.lambda(i) = 0;
                end
                if(this.freeIndex(i) && ~this.freeIndex(i+1))
                    this.lambda(i:i+2) = this.lambda(i) * this.lambdad(i:i+2);
                else
                    if (norm(this.lambda(i+1:i+2)) > this.mu * this.lambda(i))
                        scale = this.mu * this.lambda(i) / norm(this.lambda(i+1:i+2));
                        this.lambda(i+1:i+2) = scale * this.lambda(i+1:i+2);
                    end
                end
            end
        end

        %%
        function feasible = update_cg(this, alpha)
            feasible = true;
            this.lambda(this.freeIndex) = this.lambda(this.freeIndex) + alpha * this.d_cg(this.freeIndex);
            for i = 1:3:3*this.contactNum
                if(this.freeIndex(i) && this.lambda(i) < 0)
                    feasible = false;
                end
    
                if (this.freeIndex(i+1) && norm(this.lambda(i+1:i+2)) > this.mu * this.lambda(i))
                    feasible = false;
                end
            end
        end
    end

    methods(Static)
        %% Ray-Cone Intersection
        function t = rayConeIntersection(x, g, mu)
            gnorm = norm(g);
            g = g / gnorm;
        
            A = g(2)^2 + g(3)^2 - g(1)^2 * mu^2;
            B = 2 * (x(2)*g(2) + x(3)*g(3) - x(1)*g(1) * mu^2);
            C = x(2)^2 + x(3)^2 - x(1)^2 * mu^2;
        
            if(abs(C)<1e-9)
                C = 0;
            end
        
            if(abs(C)<1e-9 && gnorm < 1e-9)
                t = 0;
                return;
            end
        
            % Solve the quadratic equation
            discriminant = B^2 - 4*A*C;
        
            if(abs(A)<1e-12)
                if(g(1)>0)
                    t = Inf;
                    return;
                end
        
                if(abs(B)<1e-12)
                    t = - x(1) / g(1);
                else
                    t = -C/B;
                end
                t = t * gnorm;
                return;
            end
        
            if discriminant < 1e-12
                if g(1)>0
                    if(norm(g(2:3)) > mu*g(1))
                        t = 0;
                    else
                        t = Inf;
                    end
                else
                    t = -B / (2 * A);
                end
                t = t * gnorm;
                return;
            end
        
            % Compute the roots
            if(A>0)
                t1 = (-B - sqrt(discriminant)) / (2 * A);
                t2 = (-B + sqrt(discriminant)) / (2 * A);
            else
                t2 = (-B - sqrt(discriminant)) / (2 * A);
                t1 = (-B + sqrt(discriminant)) / (2 * A);
            end
        
            if(t1 ==0)
                if(x(1) + t2*g(1) >0)
                    t = t2;
                else
                    t = t1;
                end
                t = t * gnorm;
                return;
            end
        
            if(t2 ==0)
                if(x(1) + t1*g(1) >0)
                    t = t2;
                else
                    t = Inf;
                end
                t = t * gnorm;
                return;
            end
            
            
            if(t1>0 && t2 >0)
                t = min(t1,t2);
            elseif(t1<0 && t2<0)
                t = Inf;
            else
                t = max(t1,t2);
            end
            t = t * gnorm;
        end
    end
end