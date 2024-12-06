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
            lambda = x;
            save('GS_lambda.mat',"lambda");
        end

        function x = Chel(this, A, b, mu)
            n = size(b,1);
            x = zeros(n,1);
            y = pinv(A) * b;
            AinvA= pinv(A) * A;

            for iter = 1:this.itermax
                for i = 1 : n
                    ri = y(i) - AinvA(i,:)*x;
                    x(i) = x(i) + ri / AinvA(i,i);
                    
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
        end

        function x = Mix_GS_CG(this, A, b, mu)
            gsItermax = 25;
            cgItermax = 50;
            GSsolver = ConstraintSolver(gsItermax,this.tol);
            CGsolver = ConstraintSolver(cgItermax,this.tol);
            n = size(b,1);
            lambdas = zeros(n,1);
            iters = 1;
            freeInds = true(n,1);

            freeLambdas = CGsolver.Conjugate_Gradient(A, b, mu, lambdas, freeInds);
            lambdas(freeInds) = freeLambdas;
            this.rs(iters:iters+cgItermax - 1) = CGsolver.rs;
            iters = iters + cgItermax;
            while(iters < this.itermax)
                %lambdas = zeros(n,1);
                %{
                %freeInds = true(n,1);
                for i = 1 : n
                    if(mod(i,3) == 0)
                        if(lambdas(i-2) <= 0)
                            lambdas(i-2) = 0;
                            freeInds(i-2) = false;
                        end
                        if(norm([lambdas(i-1) lambdas(i)]) > mu * lambdas(i-2))
                            scale = mu * lambdas(i-2) / norm([lambdas(i-1) lambdas(i)]);
                            lambdas(i-1) = scale * lambdas(i-1);
                            lambdas(i) = scale * lambdas(i);
                            freeInds(i-1) = false;
                            freeInds(i) = false;
                        end
                    end
                end
                %}

                lambdas = GSsolver.Gauss_Sidiel(A, b, mu, lambdas);
                %freeInds = true(n,1);
                for i = 1 : n
                    if(mod(i,3) == 0)
                        if(lambdas(i-2) <= 0)
                            freeInds(i-2) = false;
                        end
                        if(norm([lambdas(i-1) lambdas(i)]) >= mu * lambdas(i-2))
                            freeInds(i-1) = false;
                            freeInds(i) = false;
                        end
                    end
                end
                this.rs(iters:iters+gsItermax - 1) = GSsolver.rs;
                iters = iters + gsItermax;
            end
            x = lambdas;
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
            A = A(1:3:end,1:3:end);
            Asp = Asp(1:3:end,1:3:end);
            b = b(1:3:end);
            for i = 1:length(blocks)
                blocks{i} = (1:8) + (i-1)*8;
            end

            Aspinv = Asp;
            for i = 2:5
                Aspinv(blocks{i},blocks{i})= A(blocks{i},blocks{i}) - A(blocks{i-1},blocks{i})' * pinv(Aspinv(blocks{i-1},blocks{i-1})) *  A(blocks{i-1},blocks{i});
            end
            Asp = Aspinv;

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
                        if(x(j) < 0)
                            x(j) = 0;
                        end
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