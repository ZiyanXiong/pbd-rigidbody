classdef ConstraintSolver < handle
    properties
        itermax
        rs
        tol
    end

    methods
        function this = ConstraintSolver(itermax, tol)
            this.itermax = itermax;
            this.tol = tol;
            this.rs = zeros(this.itermax,1);
        end

        function x = Gauss_Sidiel(this, A, b, mu, x0)
            n = size(b,1);
            if nargin < 5
                x = zeros(n,1);
            else
                x = x0;
            end
            for iter = 1:this.itermax
                %{
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
                %}

                db = b-A*x;
                dx = zeros(n,1);
                for normal_iter = 1:30
                    for i = 1 : 3 : n
                        ri = db(i) - A(i,:)*dx;
                        dx(i) = dx(i) + ri / A(i,i);
                        if(dx(i) < -x(i))
                            dx(i) = -x(i);
                        end
                    end
                end
                x = x+dx;

                db = b - A*x;
                dx = zeros(n,1);
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
                x = x+dx;
                this.rs(iter) = norm(b - A*(x));
            end

            lambda = x;
            save('GS_lambda.mat',"lambda");
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
            this.itermax = 1000;
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
            for iter = 1:this.itermax
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
                this.rs(iter) = norm(b - A*x);
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
            [x, f, exitflag, output, lambda]= cone_gpqp(A,-b,l,u,x,options,mu);
            %[x, f, exitflag, output]= matrixfree_gpqp(A,-b,l,u,x,options,mu);
            if(output.iterations>200)
                disp(output.iterations);
            end
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

        function x = Preconditioned_Conjugate_Gradient(this, A, b, Asp, mu)
            %Asp = diag(diag(A));
            %Asp = eye(size(A));
            [x,~,relres,iters] = pcg(A,b,this.tol,this.itermax, @(x) solvMinvx(Asp, x));
            this.rs = relres;

            function y = solvMinvx(M,x)
                Minv = pinv(M);
                y = Minv*x;
            end
        end

        function x = Shock_Propagation(this, A, b, Asp, blocks, mu)
            %{
            A = A(1:3:end,1:3:end);
            Asp = Asp(1:3:end,1:3:end);
            b = b(1:3:end);
            for i = 1:length(blocks)
                blocks{i} = (1:8) + (i-1)*8;
            end
            blocks{end} = [blocks{end}, (1:4) + i*8];
            
            Aspinv = Asp;
            for i = 2:5
                Aspinv(blocks{i},blocks{i})= A(blocks{i},blocks{i}) - A(blocks{i-1},blocks{i})' * pinv(Aspinv(blocks{i-1},blocks{i-1})) *  A(blocks{i-1},blocks{i});
            end
            Asp = Aspinv;
            %}

            n = size(b,1);
            x = zeros(n,1);

            AspT = Asp';
            for iter = 1:this.itermax
                for i = 1:length(blocks)
                %r = bsp - B*x;
                    for j = blocks{i}
                        ri = b(j) - Asp(j,:) * x;
                        x(j) = x(j) + ri ./ Asp(j,j);
                        %{
                        if(x(j) < 0)
                            x(j) = 0;
                        end
                        %}
                        %{
                        if(mod(j,3) == 0)
                            if(x(j-2) < 0)
                                x(j-2) = 0;
                            end
                            if (norm([x(j-1) x(j)]) > mu * x(j-2))
                                scale = mu * x(j-2) / norm([x(j-1) x(j)]);
                                x(j-1) = scale * x(j-1);
                                x(j) = scale * x(j);
                            end
                        end
                        %}
                    end
                end
                this.rs(iter) = norm(b - Asp*x);
            end
            %x = pinv(Asp) * b;
            rsp = b - tril(Asp - AspT) * x;
            %x = pinv(AspT)*rsp;
            x = zeros(n,1);
            for iter = 1:this.itermax
                for i = 1:length(blocks)
                    for j = blocks{i}
                        ri = rsp(j) - AspT(j,:) * x;
                        x(j) = x(j) + ri ./ AspT(j,j);
                        %{
                        if(x(j) < 0)
                            x(j) = 0;
                        end
                        %}
                        
                        if(mod(j,3) == 0)
                            if(x(j-2) < 0)
                                x(j-2) = 0;
                            end
                            if (norm([x(j-1) x(j)]) > mu * x(j-2))
                                scale = mu * x(j-2) / norm([x(j-1) x(j)]);
                                x(j-1) = scale * x(j-1);
                                x(j) = scale * x(j);
                            end
                        end
                        
                    end
                end
                this.rs(iter) = norm(b - A*x);
            end  
        end


        %%
        function draw(this, name)
            semilogy(1:size(this.rs,1), this.rs, 'DisplayName',name,'linewidth',2);
            hold on;
        end
    end
end