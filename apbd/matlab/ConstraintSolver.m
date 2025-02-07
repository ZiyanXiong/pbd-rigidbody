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

        function [lambdax] = Temporal_Gauss_Sidiel(this, A, b, d, blocks, mu, substeps)
            n = size(b,1);
            layers = length(blocks);
            lambdax = zeros(n,1);
            x = zeros(n,1);
            cpv0 = -(b + d);
            d = d * substeps;
            bsub = - (cpv0 + d);
            bsub = bsub - (cpv0+A*x);
            lambdax = lambdax + x / substeps;
            for iter = 1:substeps
                for l = 1:layers
                    nind = blocks{l};
                    nind = nind(1:3:end);
                    tind = blocks{l};
                    tind = tind(~(mod(tind,3) == 1));
                    for i = nind
                        ri = bsub(i) - A(i,:)*x;
                        x(i) = x(i) + ri / A(i,i);
                        if(x(i) < 0)
                            x(i) = 0;
                        end
                    end

                    for i = tind
                        ri = bsub(i) - A(i,:)*x;
                        x(i) = x(i) + ri / A(i,i); 
                        if (mod(i,3) == 0 && norm([x(i-1) x(i)]) > mu * x(i-2))
                            scale = mu * x(i-2) / norm([x(i-1) x(i)]);
                            x(i-1) = scale * x(i-1);
                            x(i) = scale * x(i);
                        end
                    end
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

        function x = GPQP_Warmstarted_GS(this, A, b, mu)
            x_n = this.Gradient_Projection_QP(A(1:3:end,1:3:end),b(1:3:end));
            n = size(b,1);
            x0= zeros(n,1);
            x0(1:3:n) = x_n;
            %x= x0;
            x = this.Gauss_Sidiel(A,b,mu,x0);
        end

        function x = Mix_GPQP_GS(this, A, b, mu)
            GSitermax = 1;
            rs_ = zeros(1000,1);
            n = size(b,1);
            x = zeros(n,1);
            for j = 1:1000
                db = b-A*x;
                %db = db(1:3:n) - A(1:3:end,1:3:end)*x(1:3:end);
                dx_n = this.Gradient_Projection_QP(A(1:3:end,1:3:end),db(1:3:n),-x(1:3:n));
                x(1:3:n) = x(1:3:n) + dx_n;
                %x= x0;
                db = b - A*x;
                dx = zeros(n,1);
                for iter = 1:GSitermax
                    for i = 1 : n
                        if(mod(i,3) == 1)
                            continue;
                        end
                        ri = db(i) - A(i,:)*dx;
                        dx(i) = dx(i) + ri / A(i,i);
                        
                        %Projection
                        if(mod(i,3) == 0)
                            if(x(i-2) + dx(i-2) < 0)
                                dx(i-2) = -x(i-2);
                            end
                            if (norm([x(i-1) x(i)] + [dx(i-1) dx(i)]) > mu * x(i-2))
                                scale = mu * x(i-2) / norm([x(i-1) x(i)] + [dx(i-1) dx(i)]);
                                dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                dx(i) = scale * (x(i) + dx(i)) - x(i);
                            end
                        end
                    end
                    this.rs(iter) = norm(b - A*(x+dx));
                end
                x = x + dx;
                rs_((j-1)*GSitermax+1:j*GSitermax) = this.rs(1:GSitermax);
            end
            this.rs = rs_;
        end

        function x = Staggered(this, A, b, mu)
            options.ProjectionMethod = 'none';
            options.MaxIterations = 100;
            options.Tolerance = 1e-9;
            n = length(b);
            l_n = zeros(n/3,1);
            u_n = inf(n/3,1);
            l_t = -inf(2*n/3,1);
            u_t = inf(2*n/3,1);
            x_n = zeros(n/3,1);
            x_t = zeros(2*n/3,1);
            tind = zeros(2*n/3,1);
            x = zeros(n,1);
            for j = 1:length(x_n)
                tind(j*2-1:j*2) = [(j-1)*3+1 (j-1)*3+2]+1;
            end
            iter = 1;
            this.rs = zeros(this.itermax,1);
            while(iter < this.itermax)
                db = b- A*x;
                [dx_n, f, exitflag, output, lambda]= gpqp(A(1:3:n,1:3:n),-db(1:3:n),-x_n,u_n,zeros(n/3,1), options);
                %[dx_n,~,~,~,lambdaqp] = quadprog(A(1:3:n,1:3:n),-db(1:3:n),[],[],[],[],l_n,[],[]);
                x_n = x_n + dx_n;
                x(1:3:n) = x_n;
                for j = 1:length(x_n)
                    l_t(j*2-1:j*2) = -x_n(j)*mu;
                    u_t(j*2-1:j*2) = x_n(j)*mu;
                end
                l_t = l_t - x_t;
                u_t = u_t - x_t;
                db = b - A*x;
                [dx_t, f, exitflag, output, lambda]= gpqp(A(tind,tind),-db(tind),l_t,u_t,zeros(2*n/3,1), options);
                x_t = x_t + dx_t;
                x(tind) = x_t;
                this.rs(iter:iter + sum(output.cgiterations)) = norm(b - A*x);
                iter = iter + sum(output.cgiterations) + 1;
            end
        end

        function x = Semidefinite_Programming(this, L, b, mu)
            n = length(b);
            cvx_begin sdp quiet
                variable theta;
                variable x(n);
                expression X(n,n);
                expression C(3,3,n/3);
                X = [eye(size(L,2)) L'*x; x'*L theta + 2*b'*x;];
                for i = 1:3:n
                    C(:,:,(i-1)/3+1) =[mu*x(i) 0 x(i+1); 0 mu*x(i) x(i+2); x(i+1) x(i+2) mu*x(i)];
                end
                minimize( theta );
                X >= 0;
                for i = 1:n/3
                    C(:,:,i) >= 0;
                end
            cvx_end
        end

        function x = SOCP(this, L, b, mu)
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
        end

        function x = Conjugate_Gradient(this, A, b, mu, x0, freeInds)
            n = size(b,1);
            if nargin < 5
                x = zeros(n,1);
            else
                x = x0(freeInds);
            end
            if nargin < 6
                x0 = zeros(n,1);
                x = x0;
                freeInds = true(n,1);
            end
            A0 = A;
            b0 = b;
            A = A0(freeInds, freeInds);
            b = b0(freeInds) - A0(freeInds, ~freeInds) * x0(~freeInds);
            n = size(b,1);

            r = b - A * x;                % Initial residual
            p = r;                        % Initial direction        
            % Conjugate Gradient Loop
            for iter = 1:this.itermax
                Ap = A * p;               % Matrix-vector multiplication
                alpha = (r' * r) / (p' * Ap);  % Step size
                x = x + alpha * p;        % Update solution
                r_new = r - alpha * Ap;   % Update residual
        
                % Store the norm of the residual
                x0(freeInds) = x;
                this.rs(iter) = norm(b0 -A0*x0);

                % Test if the results are valid
                isValid = true;
                for i = 1 : n
                    if(mod(i,3) == 0)
                        if(x(i-2) < 0 || norm([x(i-1) x(i)]) > mu * x(i-2))
                            %isValid = false;
                        end
                    end
                end
        
                % Check for convergence
                if norm(r_new) < this.tol || ~isValid
                    this.rs(iter:end) = this.rs(iter);
                    break;
                end
        
                % Update direction
                beta = (r_new' * r_new) / (r' * r);
                p = r_new + beta * p;
        
                % Update residual and iteration counter
                r = r_new;
            end
            % Test if the results are valid
            isValid = true;
            for i = 1 : n
                if(mod(i,3) == 0)
                    if(x(i-2) < 0 || norm([x(i-1) x(i)]) > mu * x(i-2))
                        isValid = false;
                    end
                end
            end
        end
        
        function x = Cone_GPQP(this, A, b, mu)
            options.ProjectionMethod = 'direct';
            options.MaxIterations = 1000;
            options.CGMaxIterations=200;
            options.Tolerance = 1e-8;

            l = -inf(length(b),1);
            for i = 1:3:length(b)
                l(i) = 0;
            end
            u = inf(length(b),1);
            x = zeros(length(b),1);
            n = length(b);

            [x_n, f, exitflag, output, lambda]= gpqp(A(1:3:n,1:3:n),-b(1:3:n),x(1:3:n),u(1:3:n),zeros(n/3,1), options);
            x(1:3:n) = x_n;

            [x, f, exitflag, output, lambda]= cone_gpqp(A,-b,l,u,x,options,mu);

            this.rs = zeros(this.itermax,1);
            iter = 1;
            for i = 1:length(output.cgiterations)
                this.rs(iter:iter+output.cgiterations(i)) = output.rs(i);
                iter = iter+output.cgiterations(i) + 1;
                if(iter > this.itermax)
                    break;
                end
            end

            %[x, f, exitflag, output]= matrixfree_gpqp(A,-b,l,u,x,options,mu);
            %{
            if(output.iterations>200)
                disp(output.iterations);
            end
            %}
            %[xqp,~,~,~,lambdaqp] = quadprog(A,-b,[],[],[],[],l,[],[]);
        end

        function x = Gradient_Projection_QP(this, A, b,l)
            options.ProjectionMethod = 'none';
            options.MaxIterations = 10;
            options.CGMaxIterations=100;
            options.Tolerance = 1e-9;
            if nargin < 4
                l = zeros(length(b),1);
            end
            u = inf(length(b),1);
            x = zeros(length(b),1);
            [x, f, exitflag, output, lambda]= gpqp(A,-b,l,u,x, options);
            %[xqp,~,~,~,lambdaqp] = quadprog(A,-b,[],[],[],[],l,[],[]);
        end

        function x = Preconditioned_Conjugate_Gradient(this, A, b, A_tilde)
            %Asp = diag(diag(A));
            %Asp = eye(size(A));
            [x,~,relres,iters] = pcg(A,b,this.tol,this.itermax, @(x) solvMinvx(A_tilde, x));
            this.rs = relres;

            opts.diagcomp = 1e-3;
            Li = ichol(sparse(A),opts);
            [x,~,relres,iters] = pcg(A,b,this.tol,this.itermax,Li,Li');

            function y = solvMinvx(M,x)
                Minv = pinv(M);
                y = Minv*x;
                %y = x;
            end
        end

        function x = Shock_Propagation(this, A, b, Asp, mu)
            n = size(b,1);
            x = zeros(n,1);
            this.rs = zeros(this.itermax,1);
            AspT = Asp';
            for iter = 1:1000
                for i = 1:n
                    ri = b(i) - Asp(i,:) * x;
                    x(i) = x(i) + ri ./ Asp(i,i); 
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

            rsp = b - tril(Asp - AspT) * x;
            %x = pinv(AspT)*rsp;
            %x = zeros(n,1);
            for iter = 1001:3000
                for i = 1:n
                    ri = rsp(i) - AspT(i,:) * x;
                    x(i) = x(i) + ri ./ AspT(i,i);
        
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
            %{
            for iter = 201:this.itermax
                for i = 1 : n
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
                this.rs(iter) = norm(b - A*x);
            end
            %}
        end

        function x = Shock_Propagation_Mix(this, L,Lsp, b, blocks, mu)
            n = size(b,1);
            x = zeros(n,1);
            this.rs = zeros(this.itermax*2,1);
            layers = length(blocks);
            A = L*L';
            %Lsp(73:end,:) = L(73:end,:);
            Asp = L*Lsp';
            AspT = Asp';
            for outer_iter = 1:6
                r = b - A*x;
                dx = zeros(n,1);
                for l = 1:layers
                    for iter = 1:250
                        for i = blocks{l}
                            ri = r(i) - Asp(i,:) * dx;
                            dx(i) = dx(i) + ri ./ Asp(i,i); 
                            %{
                            if(mod(i,3) == 0)
                                if(x(i-2) + dx(i-2) < 0)
                                    dx(i-2) = -x(i-2);
                                end
                                if (norm([x(i-1) + dx(i-1) x(i) + dx(i)]) > mu * (x(i-2)+dx(i-2)))
                                    scale =  mu * (x(i-2)+dx(i-2)) / norm([x(i-1) + dx(i-1) x(i) + dx(i)]);
                                    dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                    dx(i) = scale * (x(i) + dx(i)) - x(i);
                                end
                            end
                            %}
                        end
                        this.rs((outer_iter-1)*500 + iter) = norm(b - A*(x+dx));
                    end
                end
                %{
                for iter = 1:500
                    for i = 49:n
                        ri = r(i) - Asp(i,:) * dx;
                        dx(i) = dx(i) + ri ./ Asp(i,i); 
                        if(mod(i,3) == 0)
                            if(x(i-2) + dx(i-2) < 0)
                                dx(i-2) = -x(i-2);
                            end
                            if (norm([x(i-1) + dx(i-1) x(i) + dx(i)]) > mu * (x(i-2)+dx(i-2)))
                                scale =  mu * (x(i-2)+dx(i-2)) / norm([x(i-1) + dx(i-1) x(i) + dx(i)]);
                                dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                dx(i) = scale * (x(i) + dx(i)) - x(i);
                            end
                        end
                    end
                    this.rs((outer_iter-1)*1000 + iter) = norm(b - A*(x+dx));
                end
                %}
                rsp = r - tril(Asp - AspT) * dx;
                %rsp = b - (A - AspT) * x;
                %x = pinv(AspT)*rsp;
                %x = zeros(n,1);
                for l = layers:-1:1
                    for iter = 251:500
                        for i = blocks{l}
                            ri = rsp(i) - AspT(i,:) * dx;
                            dx(i) = dx(i) + ri ./ AspT(i,i);
                
                            if(mod(i,3) == 0)
                                if(x(i-2) + dx(i-2) < 0)
                                    dx(i-2) = -x(i-2);
                                end
                                if (norm([x(i-1) + dx(i-1) x(i) + dx(i)]) > mu * (x(i-2)+dx(i-2)))
                                    scale =  mu * (x(i-2)+dx(i-2)) / norm([x(i-1) + dx(i-1) x(i) + dx(i)]);
                                    dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                    dx(i) = scale * (x(i) + dx(i)) - x(i);
                                end
                            end
                        end
                        this.rs((outer_iter-1)*500 + iter) = norm(b - A*(x+dx));
                    end
                end
                %{
                for iter = 501:1000
                    for i = 1:48
                        ri = rsp(i) - AspT(i,:) * dx;
                        dx(i) = dx(i) + ri ./ AspT(i,i);
            
                        if(mod(i,3) == 0)
                            if(x(i-2) + dx(i-2) < 0)
                                dx(i-2) = -x(i-2);
                            end
                            if (norm([x(i-1) + dx(i-1) x(i) + dx(i)]) > mu * (x(i-2)+dx(i-2)))
                                scale =  mu * (x(i-2)+dx(i-2)) / norm([x(i-1) + dx(i-1) x(i) + dx(i)]);
                                dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                dx(i) = scale * (x(i) + dx(i)) - x(i);
                            end
                        end
                    end
                    this.rs((outer_iter-1)*1000 + iter) = norm(b - A*(x+dx));
                end
                %}
                x = x+ dx;
                %r = b - A*x;
                %unsatisfied_ind = abs(r)>1e-6;
                %Asp(unsatisfied_ind,unsatisfied_ind) = A(unsatisfied_ind, unsatisfied_ind);
                %Asp = A;
                %AspT = Asp';
            end
        end

        function x = Shock_Propagation_lbl(this, A, Asp, b, d, blocks, mu)
            n = size(b,1);
            x = zeros(n,1);
            upiterMax = 150;
            downiterMax = 150;
            substeps = 150;
            iterTotal = 0;
            upwardSuccess = true;
            downwardSuccess = true;
            this.rs = zeros(this.itermax,1);
            AspT = Asp';
            layers = length(blocks);
            r = b - A * x;
            dx = zeros(n,1);
            nc = 0;

            for l = 1:layers
                nind = blocks{l};
                nind = nind(1:3:end);
                tind = blocks{l};
                tind = tind(~(mod(tind,3) == 1));
                for iter = 1:upiterMax
                    rsl0 = r(blocks{l})-Asp(blocks{l},:)*dx;

                    for i = nind
                        ri = r(i) - Asp(i,:) * dx;
                        dx(i) = dx(i) + ri ./ Asp(i,i); 
                        
                        if(mod(i,3) == 1)
                            if(x(i) + dx(i) < 0)
                                dx(i) = -x(i);
                            end
                        end
                        
                        nc = nc + 1;
                        if(mod(nc,n)==0)
                            rsg = b - A*(x+dx);
                            this.rs(nc/n) = norm(rsg(rsg>0));
                        end
                    end

                    for i = tind
                        ri = r(i) - Asp(i,:) * dx;
                        dx(i) = dx(i) + ri ./ Asp(i,i); 
                        
                        if(mod(i,3) == 0)
                            if(x(i-2) + dx(i-2) < 0)
                                dx(i-2) = -x(i-2);
                            end
                            if (norm([x(i-1) + dx(i-1) x(i) + dx(i)]) > mu * (x(i-2)+dx(i-2)))
                                scale =  mu * (x(i-2)+dx(i-2)) / norm([x(i-1) + dx(i-1) x(i) + dx(i)]);
                                dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                dx(i) = scale * (x(i) + dx(i)) - x(i);
                            end
                        end
                        
                        nc = nc + 1;
                        if(mod(nc,n)==0)
                            rsg = b - A*(x+dx);
                            this.rs(nc/n) = norm(rsg(rsg>0));
                        end
                    end

                    rsl = r(blocks{l})-Asp(blocks{l},:)*dx;
                    if(norm(rsl-rsl0) < 1e-4)
                        break;
                    end
                end
                rsln = rsl(1:3:end);
                deltax = -Asp(blocks{l},:)*dx;
                if(l>1 && l < layers && ~all(rsln(deltax(1:3:end)<-2e-2) > -1e-2))
                    upwardSuccess = false;
                    break;
                end
            end
            rsp = r - tril(Asp - AspT) * dx;
            %rsp = b - (A - AspT) * x;
            %x = pinv(AspT)*rsp;
            %x = zeros(n,1);
            if(upwardSuccess)
                for l = layers:-1:1
                    nind = blocks{l};
                    nind = nind(1:3:end);
                    tind = blocks{l};
                    tind = tind(~(mod(tind,3) == 1));
                    for iter = 1:downiterMax
                        rsl0 = rsp(blocks{l})-AspT(blocks{l},:)*dx;
    
                        for i = nind
                            ri = rsp(i) - AspT(i,:) * dx;
                            dx(i) = dx(i) + ri ./ AspT(i,i); 
                            if(mod(i,3) == 1)
                                if(x(i) + dx(i) < 0)
                                    dx(i) = -x(i);
                                end
                            end
                            nc = nc + 1;
                            if(mod(nc,n)==0)
                                rsg = b - A*(x+dx);
                                this.rs(nc/n) = norm(rsg(rsg>0));
                            end
                        end
    
                        for i = tind
                            ri = rsp(i) - AspT(i,:) * dx;
                            dx(i) = dx(i) + ri ./ AspT(i,i); 
                            
                            if(mod(i,3) == 0)
                                if(x(i-2) + dx(i-2) < 0)
                                    dx(i-2) = -x(i-2);
                                end
                                if (norm([x(i-1) + dx(i-1) x(i) + dx(i)]) > mu * (x(i-2)+dx(i-2)))
                                    scale =  mu * (x(i-2)+dx(i-2)) / norm([x(i-1) + dx(i-1) x(i) + dx(i)]);
                                    dx(i-1) = scale * (x(i-1) + dx(i-1)) - x(i-1);
                                    dx(i) = scale * (x(i) + dx(i)) - x(i);
                                end
                            end
                            
                            nc = nc + 1;
                            if(mod(nc,n)==0)
                                rsg = b - A*(x+dx);
                                this.rs(nc/n) = norm(rsg(rsg>0));
                            end
                        end
                        rsl = rsp(blocks{l})-AspT(blocks{l},:)*dx;
                        if(norm(rsl-rsl0) < 1e-4)
                            break;
                        end
                    end
                    rsln = rsl(1:3:end);
                    deltax = -AspT(blocks{l},:)*dx;
                    if(l>1 && l < layers && ~all(rsln(deltax(1:3:end)<-2e-2) > -1e-2))
                        downwardSuccess = false;
                        break;
                    end
                end
            end
            x = x + dx;
            iterTotal = max(ceil(nc/n),1);
            this.itercount = iterTotal;
            rsg = b - A*x;
            this.rs(iterTotal) = norm(rsg(rsg>0));
            
            if(~downwardSuccess||~upwardSuccess)
                %{
               for iter = 1:GSiterMax
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
                    this.rs(iterTotal + iter) = norm(r(r>0));
                    if(abs(this.rs(iterTotal + iter) - this.rs(iterTotal + iter-1)) < 1e-2)
                        break;
                    end
               end 
                %}

                lambdax = zeros(n,1);
                x = zeros(n,1);
                cpv0 = -(b + d);
                d = d * substeps;
                bsub = - (cpv0 + d);
                bsub = bsub - (cpv0+A*x);
                lambdax = lambdax + x / substeps;
                for iter = 1:substeps
                    for l = 1:layers
                        nind = blocks{l};
                        nind = nind(1:3:end);
                        tind = blocks{l};
                        tind = tind(~(mod(tind,3) == 1));
                        for i = nind
                            ri = bsub(i) - A(i,:)*x;
                            x(i) = x(i) + ri / A(i,i);
                            if(x(i) < 0)
                                x(i) = 0;
                            end
                        end
    
                        for i = tind
                            ri = bsub(i) - A(i,:)*x;
                            x(i) = x(i) + ri / A(i,i); 
                            if (mod(i,3) == 0 && norm([x(i-1) x(i)]) > mu * x(i-2))
                                scale = mu * x(i-2) / norm([x(i-1) x(i)]);
                                x(i-1) = scale * x(i-1);
                                x(i) = scale * x(i);
                            end
                        end
                    end
                    bsub = bsub - (cpv0+A*x);
                    lambdax = lambdax + x / substeps;
                    rx = b - A*lambdax;
                    this.rs(iterTotal + iter) = norm(rx(rx>0));
                end
                x = lambdax;
                this.itercount = iterTotal + substeps;
            end
        end

        %%
        function draw(this, name)
            semilogy(1:size(this.rs,1), this.rs, 'DisplayName',name,'linewidth',2);
            hold on;
        end
    end
end