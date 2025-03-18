classdef ConstraintSolver < handle
    properties
        itermax
        rs
        tol
        itercount
    end

    methods
        function this = ConstraintSolver(itermax, tol)
            this.itermax = itermax;
            this.tol = tol;
            this.rs = zeros(this.itermax,1);
            this.itercount = 0;
        end

        function x = Gauss_Sidiel(this, A, b, mu, x0)
            n = size(b,1);
            if nargin < 5
                x = zeros(n,1);
            else
                x = x0;
            end
            this.rs = zeros(this.itermax,1);
            for iter = 1:this.itermax
                
                for i = 1 : 3 : n
                    ri = b(i) - A(i,:)*x;
                    x(i) = x(i) + ri / A(i,i);
                    if(x(i) < 0)
                        x(i) = 0;
                    end
                end
                
                for i = 1 : n
                    if(mod(i,3) == 1)
                        continue;
                    end

                    ri = b(i) - A(i,:)*x;
                    x(i) = x(i) + ri / A(i,i);
                   
                    %Projection
                    if(mod(i,3) == 0)
                        if(x(i-2) < 0)
                            x(i-2) = 0;
                        end
                        if (norm([x(i-1) x(i)]) > mu * x(i-2))
                            scale = mu * x(i-2) / norm([x(i-1) x(i)]);
                            x(i-1) = scale * x(i-1);
                            x(i) = scale * x(i);
                        end
                    end
                   
                end
                r = b - A*x;
                this.rs(iter) = norm(r(r>0));
            end
            this.itercount = iter;
            %lambda = x;
            %save('GS_lambda.mat',"lambda");
        end

        function [lambdax, x] = Temporal_Gauss_Sidiel(this, A, b, d, contactConstraintEndInd, mu, substeps, x0)
            n = size(b,1);
            lambdax = zeros(n,1);
            if(nargin <8)
                x = zeros(n,1);
            else
                x = x0;
            end
            cpv0 = -(b + d);
            d = d * substeps;
            bsub = - (cpv0 + d);
            for iter = 1:substeps
                nind = 1:3:contactConstraintEndInd;
                tind = 2:3:contactConstraintEndInd;
                for i = nind
                    ri = bsub(i) - A(i,:)*x;
                    x(i) = x(i) + ri / A(i,i);
                    if(x(i) < 0)
                        x(i) = 0;
                    end
                end

                for i = tind
                    ri = bsub(i:i+1) - A(i:i+1,:) * x;
                    x(i) = x(i) + ri(1) ./ A(i,i); 
                    x(i+1) = x(i+1) + ri(2) ./ A(i+1,i+1);

                    if (norm([x(i) x(i+1)]) > mu * x(i-1))
                        scale =  mu * (x(i-1)) / norm([x(i) x(i+1)]);
                        x(i) = scale * x(i);
                        x(i+1) = scale * x(i+1);
                    end
                end

                for i = contactConstraintEndInd + 1 : n
                    ri = bsub(i) - A(i,:)*x;
                    x(i) = x(i) + ri / A(i,i);
                end
                bsub = bsub - (cpv0+A*x);
                lambdax = lambdax + x / substeps;
                rx = b - A*lambdax;
                this.rs(iter) = norm(rx(rx>0));
            end
            this.itercount = substeps;
            %lambda = x;
            %save('GS_lambda.mat',"lambda");
        end

        function [x, xv] = SOCP(this, L, b, mu)
            n = length(b);
            dsc = zeros(n+1,1);
            dsc(end) =1;
            gamma = -1;
            Asc = [L' zeros(size(L,2),1); zeros(1,size(L,1)) 1];
            fsc = [-b; 1];
            cvx_begin quiet
                variable t;
                variable x(n);
                expression u(n+1);
                u = [x; t];
                minimize( fsc'*u );
                norm( Asc*u ) <= dsc'*u-gamma;
                for i = 1 : 3 : n
                    norm( u(i+1:i+2) ) <= mu*u(i);
                end
            cvx_end
            xv = x;
            this.itercount = cvx_slvitr;
            this.rs = zeros(this.itermax,1);
        end


        function [x, xv] = Cone_GPQP(this, A, b, mu, contactConstraintEndInd)
            options.ProjectionMethod = 'direct';
            options.MaxIterations = 1000;
            options.CGMaxIterations=200;
            options.Tolerance = 1e-6;

            l = -inf(length(b),1);
            u = inf(length(b),1);
            for i = 1:3:contactConstraintEndInd
                l(i) = 0;
            end
            x = zeros(length(b),1);

            [x, f, exitflag, output, lambda]= cone_gpqp(A,-b,l,u,x,options,mu,contactConstraintEndInd);
            xv = x;
            this.rs = zeros(this.itermax,1);
            iter = 1;
            for i = 1:output.iterations
                this.rs(iter:iter+output.cgiterations(i)-1) = output.rs(i);
                iter = iter+output.cgiterations(i);
                if(iter > this.itermax)
                    break;
                end
            end
            this.itercount = iter-1;
        end

        function [x, xv] = Staggered(this, A, b, mu)
            options.ProjectionMethod = 'none';
            options.MaxIterations = 100;
            options.Tolerance = 1e-9;
            n = length(b);
            l = zeros(n,1);
            u = inf(n,1);
            x = zeros(n,1);
            for i = 1:3:n
                l(i+1:i+2) = -x(i)*mu;
                u(i+1:i+2) = x(i)*mu;
            end
            nind = false(n,1);
            nind(1:3:n) = true;
            tind = ~nind;
            iter = 1;
            this.rs = zeros(this.itermax,1);
            %{
            while(iter < this.itermax)
                bn = b(nind) - A(nind,tind)*x(tind);
                [xn, f, exitflag, output, lambda]= gpqp(A(nind,nind),-bn,l(nind),u(nind),x(nind), options);
                %[dx_n,~,~,~,lambdaqp] = quadprog(A(1:3:n,1:3:n),-db(1:3:n),[],[],[],[],l_n,[],[]);
                x(nind) = xn;
                for i = 1:3:n
                    l(i+1:i+2) = -x(i)*mu;
                    u(i+1:i+2) = x(i)*mu;
                end

                bt = b(tind) - A(tind,nind)*x(nind);
                [xt, f, exitflag, output, lambda]= gpqp(A(tind,tind),-bt,l(tind),u(tind),x(tind), options);
                x(tind) = xt;

                this.rs(iter) = norm(b - A*x);
                iter = iter + 1;
            end
            %}
            
            [x, f, exitflag, output, lambda]= gpqp_staggered(A,-b,l,u,x,options,mu);
            this.rs = zeros(this.itermax,1);
            iter = 1;
            for i = 1:output.iterations
                this.rs(iter:iter+output.cgiterations(i)-1) = output.rs(i);
                iter = iter+output.cgiterations(i);
                if(iter > this.itermax)
                    break;
                end
            end
            
            xv = x;
            this.itercount = iter-1;
        end

        %%
        function draw(this, name)
            semilogy(1:size(this.rs,1), this.rs, 'DisplayName',name,'linewidth',2);
            hold on;
        end
    end
end