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

        function [lambdax, x] = Temporal_Gauss_Sidiel(this, A, b, d, blocks, mu, substeps, x0)
            n = size(b,1);
            layers = length(blocks);
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
                for l = 1:layers
                    for np = 1:length(blocks{l})
                        p = blocks{l}{np};
                        for i = p(1:3:end)
                            ri = bsub(i) - A(i,:) * x;
                            x(i) = x(i) + ri ./ A(i,i);
                            if(x(i) < 0)
                                x(i) = 0;
                            end
                        end
                    end

                    for np = 1:length(blocks{l})
                        p = blocks{l}{np};
                        for i = p(2:3:end)
                            ri = bsub(i:i+1) - A(i:i+1,:) * x;
                            x(i) = x(i) + ri(1) ./ A(i,i); 
                            x(i+1) = x(i+1) + ri(2) ./ A(i+1,i+1);
    
                            if (norm([x(i) x(i+1)]) > mu * x(i-1))
                                scale =  mu * (x(i-1)) / norm([x(i) x(i+1)]);
                                x(i) = scale * x(i);
                                x(i+1) = scale * x(i+1);
                            end
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

        function [lambdax, x] = Gauss_Sidiel_SP(this, A, Asp, b, blocks, mu)
            n = size(b,1);
            layers = length(blocks);
            %lambdax = zeros(n,1);
            x = zeros(n,1);
            Asp= Asp+(triu(Asp' - Asp));
            AspD = A - Asp;

            for iter = 1:1
                bsp = b - AspD*x;
                for l = 1:layers
                    for layer_iter = 1 : 75
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(1:3:end)
                                ri = bsp(i) - Asp(i,:) * x;
                                x(i) = x(i) + ri ./ Asp(i,i);
                                if(x(i) < 0)
                                    x(i) = 0;
                                end
                            end
                        end
    
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(2:3:end)
                                ri = bsp(i:i+1) - Asp(i:i+1,:) * x;
                                x(i) = x(i) + ri(1) ./ Asp(i,i); 
                                x(i+1) = x(i+1) + ri(2) ./ Asp(i+1,i+1);
        
                                if (norm([x(i) x(i+1)]) > mu * x(i-1))
                                    scale =  mu * (x(i-1)) / norm([x(i) x(i+1)]);
                                    x(i) = scale * x(i);
                                    x(i+1) = scale * x(i+1);
                                end
                            end
                        end
                    end
                end
                
                for l = layers:-1:1
                    for layer_iter = 1 : 75
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(1:3:end)
                                ri = bsp(i) - Asp(i,:) * x;
                                x(i) = x(i) + ri ./ Asp(i,i);
                                if(x(i) < 0)
                                    x(i) = 0;
                                end
                            end
                        end
    
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(2:3:end)
                                ri = bsp(i:i+1) - Asp(i:i+1,:) * x;
                                x(i) = x(i) + ri(1) ./ Asp(i,i); 
                                x(i+1) = x(i+1) + ri(2) ./ Asp(i+1,i+1);
        
                                if (norm([x(i) x(i+1)]) > mu * x(i-1))
                                    scale =  mu * (x(i-1)) / norm([x(i) x(i+1)]);
                                    x(i) = scale * x(i);
                                    x(i+1) = scale * x(i+1);
                                end
                            end
                        end
                    end
                end

                rx = b - A*x;
                this.rs(iter) = norm(rx(rx>0));
            end

            this.itercount = this.itermax;
            lambdax = x;
            %save('GS_lambda.mat',"lambda");
        end

        function [lambdax] = Temporal_Gauss_Sidiel_Joints(this, A, b, substeps)
            n = size(b,1);
            lambdax = zeros(n,1);
            x = zeros(n,1);
            d = 0;
            cpv0 = -(b + d);
            d = d * substeps;
            bsub = - (cpv0 + d);
            bsub = bsub - (cpv0+A*x);
            lambdax = lambdax + x / substeps;
            for iter = 1:substeps
                for i = 1:n
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

        function [x, lambdav] = Shock_Propagation_lbl(this, A, Asp, b, d, blocks, mu)
            n = size(b,1);
            x = zeros(n,1);
            upiterMax = 75;
            downiterMax = 75;
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
                blockinds = [blocks{l}{:}];
                for iter = 1:upiterMax
                    rsl0 = r(blockinds)-Asp(blockinds,:)*dx;
                    for np = 1:length(blocks{l})
                        p = blocks{l}{np};
                        
                        rn = r(p(1:3:end)) - Asp(p(1:3:end),:) * dx;
                        dxn = rn ./ diag(Asp(p(1:3:end),p(1:3:end)));
                        dxn((rn<=0 & dx(p(1:3:end))==0)) = 0;
                        dxn = dxn*0.5;
                        %Line search for Jacobi
                        lsIter = 1;
                        while(norm(rn - Asp(p(1:3:end),p(1:3:end))*dxn) > norm(rn))
                            dxn = dxn*0.5;
                            lsIter = lsIter + 1;
                            if(lsIter > 10)
                                dxn = dxn*0;
                                break;
                            end
                        end

                        
                        dx(p(1:3:end)) = dx(p(1:3:end)) + dxn;
                        
                        for i = p(1:3:end)
                            ri = r(i) - Asp(i,:) * dx;
                            dx(i) = dx(i) + ri / Asp(i,i);
                            if(x(i) + dx(i) < 0)
                                dx(i) = -x(i);
                            end

                            nc = nc + 1;
                            if(mod(nc,n)==0)
                                rsg = b - A*(x+dx);
                                this.rs(nc/n) = norm(rsg(rsg>0));
                            end
                        end
                    end
                    
                    for np = 1:length(blocks{l})
                        p = blocks{l}{np};
                        for i = p(2:3:end)
                            ri = r(i:i+1) - Asp(i:i+1,:) * dx;
                            %dx(i:i+1) = dx(i:i+1) + Asp(i:i+1,i:i+1)\ri; 
                            dx(i) = dx(i) + ri(1) / Asp(i,i);
                            dx(i+1) = dx(i+1) + ri(2) / Asp(i+1,i+1);
                            
                            if (norm([x(i) + dx(i) x(i+1) + dx(i+1)]) > mu * (x(i-1)+dx(i-1)))
                                scale =  mu * (x(i-1)+dx(i-1)) / norm([x(i) + dx(i) x(i+1) + dx(i+1)]);
                                dx(i) = scale * (x(i) + dx(i)) - x(i);
                                dx(i+1) = scale * (x(i+1) + dx(i+1)) - x(i+1);
                            end
                            
                            nc = nc + 2;
                            if(mod(nc,n)==0)
                                rsg = b - A*(x+dx);
                                this.rs(nc/n) = norm(rsg(rsg>0));
                            end
                        end
                    end
                    
                    rsl = r(blockinds)-Asp(blockinds,:)*dx;
                    if(norm(rsl-rsl0) < 1e-6)
                        break;
                    end
                end
                rsln = rsl(1:3:end);
                rsln0 = r(blockinds(1:3:end));
                %deltax = -Asp(blockinds,:)*dx;
                %if(~all(rsln(deltax(1:3:end)<-1) > -1))
                if(any(rsln<-3) && any(rsln0(rsln<-3) > 1))
                    upwardSuccess = false;
                    break;
                end
            end
            rsp = r - tril(Asp - AspT) * dx;
            %rsp = b - (A - AspT) * x;
            %x = pinv(AspT)*rsp;
            %dx = zeros(n,1);
            if(upwardSuccess)
                for l = layers:-1:1
                    blockinds = [blocks{l}{:}];
                    for iter = 1:downiterMax
                        rsl0 = rsp(blockinds)-AspT(blockinds,:)*dx; 
    
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            
                            rn = rsp(p(1:3:end)) - AspT(p(1:3:end),:) * dx;
                            dxn = rn ./ diag(AspT(p(1:3:end),p(1:3:end)));
                            dxn((rn<=0 & dx(p(1:3:end))==0)) = 0;
                            dxn = dxn*0.5;
                            %Line search for Jacobi
                            lsIter = 1;
                            while(norm(rn - Asp(p(1:3:end),p(1:3:end))*dxn) > norm(rn))
                                dxn = dxn*0.5;
                                lsIter = lsIter + 1;
                                if(lsIter > 10)
                                    dxn = dxn*0;
                                    break;
                                end
                            end                   
                            
                            dx(p(1:3:end)) = dx(p(1:3:end)) + dxn;
                            
                            for i = p(1:3:end)                                
                                %ri = rsp(i) - AspT(i,:) * dx;
                                %dx(i) = dx(i) + ri / AspT(i,i);
                                if(x(i) + dx(i) < 0)
                                    dx(i) = -x(i);
                                end

                                nc = nc + 1;
                                if(mod(nc,n)==0)
                                    rsg = b - A*(x+dx);
                                    this.rs(nc/n) = norm(rsg(rsg>0));
                                end
                            end
                        end

                        
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(2:3:end)
                                ri = rsp(i:i+1) - AspT(i:i+1,:) * dx;
                                %dx(i:i+1) = dx(i:i+1) + Asp(i:i+1,i:i+1)\ri; 
                                dx(i) = dx(i) + ri(1) / AspT(i,i);
                                dx(i+1) = dx(i+1) + ri(2) / AspT(i+1,i+1);
        
                                if (norm([x(i) + dx(i) x(i+1) + dx(i+1)]) > mu * (x(i-1)+dx(i-1)))
                                    scale =  mu * (x(i-1)+dx(i-1)) / norm([x(i) + dx(i) x(i+1) + dx(i+1)]);
                                    dx(i) = scale * (x(i) + dx(i)) - x(i);
                                    dx(i+1) = scale * (x(i+1) + dx(i+1)) - x(i+1);
                                end
        
                                nc = nc + 2;
                                if(mod(nc,n)==0)
                                    rsg = b - A*(x+dx);
                                    this.rs(nc/n) = norm(rsg(rsg>0));
                                end
                            end
                        end
                        
                        rsl = rsp(blockinds)-AspT(blockinds,:)*dx;
                        if(norm(rsl-rsl0) < 1e-6)
                            break;
                        end
                    end
                    rsln = rsl(1:3:end);
                    rsln0 = rsp(blockinds(1:3:end));
                    %deltax = -AspT(blockinds,:)*dx;
                    %if(~all(rsln(deltax(1:3:end)<-1) > -1))
                    if(any(rsln<-3) && any(rsln0(rsln<-3) > 1))
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
            lambdav = x;
            
            if(~downwardSuccess||~upwardSuccess)
                lambdax = zeros(n,1);
                x = zeros(n,1);
                cpv0 = -(b + d);
                d = d * substeps;
                bsub = - (cpv0 + d);
                for iter = 1:substeps
                    for l = 1:layers
                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(1:3:end)
                                ri = bsub(i) - A(i,:) * x;
                                x(i) = x(i) + ri ./ A(i,i);
                                if(x(i) < 0)
                                    x(i) = 0;
                                end
                            end
                        end

                        for np = 1:length(blocks{l})
                            p = blocks{l}{np};
                            for i = p(2:3:end)
                                ri = bsub(i:i+1) - A(i:i+1,:) * x;
                                x(i) = x(i) + ri(1) ./ A(i,i); 
                                x(i+1) = x(i+1) + ri(2) ./ A(i+1,i+1);
        
                                if (norm([x(i) x(i+1)]) > mu * x(i-1))
                                    scale =  mu * (x(i-1)) / norm([x(i) x(i+1)]);
                                    x(i) = scale * x(i);
                                    x(i+1) = scale * x(i+1);
                                end
                            end
                        end
                    end
                    bsub = bsub - (cpv0+A*x);
                    lambdax = lambdax + x / substeps;
                    rx = b - A*lambdax;
                    this.rs(iterTotal + iter) = norm(rx(rx>0));
                end
                lambdav = x;
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