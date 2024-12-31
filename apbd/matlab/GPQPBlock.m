classdef GPQPBlock < handle
    properties
        ind %indices of corresponding contact
        neighbors %indices of adjacent contacts
        neighborMatrices
        blockList
        matrix
        b
        x
        xc
        l
        u
        t_bar
        g
        p
        Ax
        Gp
        activeInd

        %% variables for CG
        r_cg
        b_cg
        Minv_cg
        g_cg
        d_cg
    end

    methods
        %%
        function this = GPQPBlock(A, b, x0, ind)
            n = length(b);
            this.blockList = {};
            this.ind = ind;
            this.matrix = A(ind:ind+2,ind:ind+2);
            this.b = b(ind:ind+2);
            this.x = x0(ind:ind+2);
            this.neighbors = [];
            this.l = [0 0 0]';
            this.u = [Inf 0 0]';
            this.t_bar = zeros(3,1);
            this.g = zeros(3,1);
            this.Ax = zeros(3,1);
            this.Gp = zeros(3,1);
            this.xc = zeros(3,1);
            this.p = zeros(3,1);
            this.activeInd = false(3,1);
            this.neighborMatrices = {};
            for i=1:3:n
                if(A(ind,i)~=0 && i ~= ind)
                    this.neighbors(end+1) = (i+2)/3;
                    this.neighborMatrices(end+1) = {A(ind:ind+2,i:i+2)};
                end
            end
        end

        %%
        function compute_Ax_g(this)
            this.Ax = this.matrix * this.x;
            for i = 1:length(this.neighbors)
                this.Ax = this.Ax + this.neighborMatrices{i} * this.blockList{this.neighbors(i)}.x;
            end
            this.g = this.Ax + this.b;
        end

        %%
        function tbar_i = compute_tbar(this)
            if((abs(this.x(1)-this.u(1))<1e-12 && abs(this.g(1)) < 1e-12)  || (abs(this.g(1)) < 1e-12 && abs(this.x(1)-this.l(1))<1e-12))
                this.t_bar(1) = 0;
            elseif((this.g(1)>0))
                this.t_bar(1) = (this.x(1) - this.l(1)) / this.g(1);
            else
                this.t_bar(1) = Inf;
            end
            xi = this.x(2:3);
            gi = this.g(2:3);
            r = this.u(2);
            if(abs(norm(xi) - r)<1e-12 && norm(gi)<1e-12)
                ti = 0;
            else
                ti = (xi'*gi + sqrt((xi'*gi)^2 -gi'*gi*(xi'*xi - r^2))) / (gi'*gi);
            end
            this.t_bar(2:3) = ti;

            if(this.t_bar(1) == 0)
                this.t_bar(2:3) = 0;
            end

            tbar_i = this.t_bar;
        end

        %%
        function compute_xc(this, t)
            for i = 1:3
                if(t < this.t_bar(i))
                    this.xc(i) = this.x(i) - t * this.g(i);
                    this.activeInd(i) = false;
                else
                    this.xc(i) = this.x(i) - this.t_bar(i) * this.g(i);
                    this.activeInd(i) = true;
                end
            end
        end

        %%
        function compute_p(this,t)
            for i = 1:3
                if(t < this.t_bar(i))
                    this.p(i) = - this.g(i);
                else
                    this.p(i) = 0;
                end
            end
        end

        %%
        function compute_Gp(this)
            this.Gp = this.matrix * this.p;
            for i = 1:length(this.neighbors)
                this.Gp = this.Gp + this.neighborMatrices{i} * this.blockList{this.neighbors(i)}.p;
            end
        end

        %%
        function project(this, mu)
            if(this.x(1) < 0)
                this.x(1) = 0;
            end

            if (norm(this.x(2:3)) > mu * this.x(1))
                scale = mu * this.x(1) / norm(this.x(2:3));
                this.x(2) = scale * this.x(2);
                this.x(3) = scale * this.x(3);
            end
            this.l(2:3) = -mu*this.x(1);
            this.u(2:3) = mu*this.x(1);
        end

        %%
        function compute_b_cg(this)
            notfreeInd = this.activeInd;
            this.b_cg = this.b + this.matrix(:,notfreeInd) * this.x(notfreeInd);
            for i = 1:length(this.neighbors)
                notfreeInd = this.blockList{this.neighbors(i)}.activeInd;
                this.b_cg = this.b_cg + this.neighborMatrices{i}(:, notfreeInd) * this.blockList{this.neighbors(i)}.x(notfreeInd);
            end
        end

        %%
        function init_cg(this)
            freeInd = ~this.activeInd;
            this.Ax = this.matrix(:,freeInd) * this.x(freeInd);
            for i = 1:length(this.neighbors)
                freeInd = ~this.blockList{this.neighbors(i)}.activeInd;
                this.Ax = this.Ax + this.neighborMatrices{i}(:, freeInd) * this.blockList{this.neighbors(i)}.x(freeInd);
            end
            this.r_cg = this.Ax + this.b_cg;
            this.g_cg = this.Minv_cg .* this.r_cg;
            this.d_cg = - this.g_cg;
        end

        %%
        function compute_Ad_cg(this)
            freeInd = ~this.activeInd;
            this.Ax = this.matrix(:,freeInd) * this.d_cg(freeInd);
            for i = 1:length(this.neighbors)
                freeInd = ~this.blockList{this.neighbors(i)}.activeInd;
                this.Ax = this.Ax + this.neighborMatrices{i}(:, freeInd) * this.blockList{this.neighbors(i)}.d_cg(freeInd);
            end
        end

        %%
        function feasible = update_cg(this, alpha)
            feasible = true;
            this.x(~this.activeInd) = this.x(~this.activeInd) + alpha * this.d_cg(~this.activeInd);
            if(~this.activeInd(1) && this.x(1) < this.l(1))
                feasible = false;
            end

            if (~this.activeInd(2) && norm(this.x(2:3)) > this.u(2))
                feasible = false;
            end
        end

    end
end