%% Numerical Optimization (Jorge Nocedal and Stephen J. Wright), Springer, 2006.
%% Gradient Projection Method for QP (Algorithm 16.5)
function [x, f, exitflag, output, lambda]= cone_gpqp(G,c,l,u,x,opts, mu)
if nargin < 2
	error('At least two arguments required');
end
n = length(c);
if nargin < 3
	l = -inf(n,1);
end
if nargin < 4
	u = inf(n,1);
end
if nargin < 5
	x = zeros(n,1);
end
if nargin < 6
	opts = [];
end
if(isempty(x))
	x = zeros(n,1);
end
if(isempty(opts))
    opts.ProjectionMethod = 'direct';
    opts.MaxIterations = 100;
    opts.Tolerance = 1e-9;
end
if(isfield(opts,'ProjectionMethod'))
    projectionOpts = opts.ProjectionMethod;
else
    projectionOpts = 'direct';
end
if(isfield(opts,'MaxIterations'))
    iterNum = opts.MaxIterations;
else
    iterNum = 100;
end
if(isfield(opts,'CGMaxIterations'))
    CGiterNum = opts.CGMaxIterations;
else
    CGiterNum = 100;
end
if(isfield(opts,'Tolerance'))
    eps = opts.Tolerance;
else
    eps = 1e-9;
end

%{
H = abs(diag(G));
[x,~,~,~,~] = pcg_rs(G, -c, 1e-8, CGiterNum,H,[],x,[],[],projectionOpts,1:length(c), mu);
for i = 1:3:n
    if(x(i) < 0)
        x(i) = 0;
    end

    if (norm([x(i+1) x(i+2)]) > mu * x(i))
        scale = (mu * x(i)) / norm([x(i+1) x(i+2)]);
        x(i+1) = scale * x(i+1);
        x(i+2) = scale * x(i+2);
    end
    
end
%}

fPrev = 0.5 * x' * G * x + x' * c;
gPrev = G * x + c;
CGiterVec = [];
rs = zeros(iterNum,1);
normalIndex = true(n,1);
tangentIndex = false(n,1);
for i = 1:3:n
    normalIndex(i+1:i+2) = false;
    tangentIndex(i+1:i+2) = true;
end

for iter = 1:iterNum
    if(iter == 10)
        disp("stop");
    end
    for i = 1:3:n
        l(i+1) = -mu*x(i);
        u(i+1) = mu*x(i);
        l(i+2) = -mu*x(i);
        u(i+2) = mu*x(i);
    end
        
    [xc,xd,neqIndex] = computeCauchyPointCone(G,c,l,u,x, mu);
    x = xc;
    selectM = zeros(n,length(find(neqIndex)));
    selecti = 1;
    for i = 1:3:n
        if(neqIndex(i)&&~neqIndex(i+1))
            selectM(i:i+2, selecti) = xd(i:i+2);
            selecti = selecti + 1;
        elseif(neqIndex(i)&&neqIndex(i+1))
            selectM(i:i+2,selecti:selecti+2) = eye(3);
            selecti = selecti + 3;
        end
    end
    
    cZ = c;
    
    %H = ones(n,1);
    %H = diag(abs(diag(G)));
    %W_zz = Z' * H * Z;
    %W_zz = H(neqIndex);

    H = abs(diag(selectM'*G*selectM));
    W_zz = H;

    %xZ = pcg_rs(Z' * G * Z, -cZ,eps,CGiterNum,W_zz,[],xZ,lNeq,uNeq,projectionOpts);
    %[xZ,~,~,CGiter,~] = pcg_rs(G(neqIndex, neqIndex), -cZ,1e-9,CGiterNum,W_zz,[],xZ,lNeq,uNeq,projectionOpts,find(neqIndex));
    %[xZ,~,~,CGiter,~] = pcr1(G(neqIndex, neqIndex), -cZ, 1e-8, CGiterNum,diag(W_zz),[],xZ,[],[],projectionOpts,find(neqIndex), mu);
    [xZ,~,~,CGiter,~] = pcr1(selectM'*G*selectM, -selectM'*cZ, 1e-7, CGiterNum,diag(W_zz),[],selectM' * x,[],[],projectionOpts,find(neqIndex), mu);
    %[xZ,~,~,CGiter,~] = pcg_rs(selectM'*G*selectM, -selectM'*cZ, 1e-8, CGiterNum,W_zz,[],selectM' * x,[],[],projectionOpts,find(neqIndex), mu);
    %[xZ,~,~,CGiter,~] = minres1(selectM'*G*selectM, -selectM'*cZ, 1e-7, CGiterNum,diag(W_zz),[],selectM' * x,[],[],projectionOpts,find(neqIndex), mu);

    %R = chol(selectM'*G*selectM + 1e-6*eye(selecti - 1));
    %[xZ,~,~,CGiter,~] = minres(selectM'*G*selectM,-selectM'*cZ,1e-4,CGiterNum,R',R, selectM' * x);
    %[xZ,~,~,CGiter,~] = minres(selectM'*G*selectM + 1e-6*eye(selecti - 1),-selectM'*cZ,1e-8,CGiterNum);

    %[xZ,~,~,CGiter,~] = gs(G(neqIndex, neqIndex), -cZ, 1e-8, 20,diag(W_zz),[],xZ,lNeq,uNeq,projectionOpts,find(neqIndex), mu);
    CGiterVec = [CGiterVec, CGiter];
    
    selecti = 1;
    for i = 1:3:n
        if(neqIndex(i)&&~neqIndex(i+1))
            if(xZ(selecti)<0)
                xZ(selecti) = 0;
            end
            x(i:i+2) = xZ(selecti)*xd(i:i+2);
            selecti = selecti + 1;
        elseif(neqIndex(i)&&neqIndex(i+1))
            x(i:i+2) = xZ(selecti:selecti+2);
            selecti = selecti + 3;
        end
    end

    for i = 1:3:n
        if(x(i) < 0)
            x(i) = 0;
        end

        if (norm([x(i+1) x(i+2)]) > mu * x(i))
            scale = (mu * x(i)) / norm([x(i+1) x(i+2)]);
            x(i+1) = scale * x(i+1);
            x(i+2) = scale * x(i+2);
        end
        
    end
    
    % if satisfies the KKT conditions
    f = 0.5 * x' * G * x + x' * c;
    g = G * x + c;

    locked = true;
    for i = 1:3:n
        if(neqIndex(i)&&neqIndex(i+1))
            locked = false;
        end
    end
    if(locked)
        for i = 1:3:n
            if(neqIndex(i)&&~neqIndex(i+1))
                if (norm([g(i+1) g(i+2)]) > 1e-6 && norm([x(i+1) x(i+2)]) > 1e-6)
                    gdi = -g(i+1:i+2) / norm(g(i+1:i+2));
                    xdi = x(i+1:i+2) / norm(x(i+1:i+2));
                else
                    continue;
                end
                
                if( gdi'*xdi > 0 && gdi'*xdi < 0.98)
                    %x(i+1:i+2) = abs(gdi'*xdi) * x(i+1:i+2);
                end
            end
        end
    end
    
    rs(iter) = norm(g);
	if norm(g) < eps && sum(neqIndex) == n
		break;
	end
    if norm(g - gPrev) < eps
        break;
    end
    if norm(f-fPrev) < eps
        break;
    end
    fPrev = f;
    gPrev = g;
end

lambdaLower = G * x + c;
lambdaUpper = lambdaLower;
lambdaLower(lambdaLower < 0) = 0;
lambdaUpper(lambdaUpper > 0) = 0;
lambda.lower = lambdaLower;
lambda.upper = -lambdaUpper;

if(iter ~= iterNum)
    exitflag = 1;
else
    exitflag = 0;
end

output.iterations = iter;
output.projectionMethod = projectionOpts;
output.cgiterations = CGiterVec;
output.rs = rs;
end

%% Preconditioned CG for Reduced Systems (Algorithm 16.1 )
function [xZ, flag, relres, j, resvec] = pcg_rs(A,b,tol,maxit,M1,~,x0,lReduced,uReduced,projectionOpts,mindices,mu)
%M1 is assumed to be a diagnal matrix, so it is represented by a vector
%inside this function.
n = size(b,1);
if isempty(tol)
    tol = 1e-6;
end
if isempty(maxit)
    maxit = 100;
end
if isempty(M1)
    M1 = ones(n);
end
if isempty(projectionOpts)
    projectionOpts = 'direct';
end

M1inv =  1 ./ M1; 
xZ = x0;
rZ = A * xZ - b;
gZ = M1inv .* rZ;
dZ = -gZ;

resvec = [];
for j = 1: maxit
    if(rZ' * diag(M1) * rZ < tol)
        flag = 0;
        break;
    end
    alpha = (rZ' * gZ) / (dZ' * A * dZ);
    xZ_last = xZ;
    xZ = xZ + alpha * dZ;
    
    feasible = true;
    for i = 1:length(xZ)
        if(mod(mindices(i),3)==1 && xZ(i) < -1e-6)
            feasible = false;
        end

        if(mod(mindices(i),3)==2 && norm([xZ(i) xZ(i+1)]) > mu*xZ(i-1) + 1e-6)
            feasible = false;
        end
    end

    if(~feasible)
        flag = 1;
        break;
    end

    %{
    % terminate as soon as a bound l ≤ x ≤ u is encountered
    if(sum(xZ < lReduced) > 0  || sum(xZ > uReduced) > 0)
        switch projectionOpts
            case 'none'
                flag = 1;
                break;
            case 'alpha'
                alphal = abs((lReduced - xZ_last) ./ dZ);
                alphau = abs((uReduced - xZ_last) ./ dZ);
                alphaList = [alphal(xZ < lReduced);alphau(xZ > uReduced)];
                alpha = sign(alpha) * min(alphaList);
                xZ = xZ_last + alpha * dZ;
                flag = 1;
                break;
            case 'direct'
                xZ(xZ < lReduced) = lReduced(xZ < lReduced);
                xZ(xZ > uReduced) = uReduced(xZ > uReduced);
                flag = 1;
                break;
            otherwise
                msg = "Invalid Projection Method. Choose from 'none', 'alpha' and 'direct'";
                error(msg)
        end
    end
    %}

    rZ_plus = rZ + alpha * A * dZ;
    gZ_plus = M1inv .* rZ_plus;
    beta = (rZ_plus' * gZ_plus) / (rZ' * gZ);
    dZ = - gZ_plus + beta * dZ;
    gZ = gZ_plus;
    rZ = rZ_plus;
    resvec =[resvec, norm(b - A* xZ)];
end
%xZ(xZ < lReduced) = lReduced(xZ < lReduced);
%xZ(xZ > uReduced) = uReduced(xZ > uReduced);
relres = norm(A * xZ - b) / norm(b);
if(j == maxit)
    flag = 2;
end
end

%% Preconditioned Conjugate Residual (GNU Octave)
function [x,flag,relres,iter,resvec] = pcr1(A,b,tol,maxit,M,~,x0,~,uReduced,~,mindices, mu)
    flag = 1;
    relres = 0;
    if nargout >= 5
	    resvec = [];
    end
    
    x = x0;
    r = b - A * x;
    p = M \ r;
    b_bot_old = 1;
    s_old = zeros(size(x));
    p_old = s_old;
    q_old = p_old;
    q = A * p;
    
    normb = norm(b);
    for iter = 0 : maxit
	    % Convergence check
	    normr0 = norm(r);
	    if nargout >= 5
		    resvec(end+1) = normr0; %#ok<AGROW>
	    end
	    if normr0 < tol*normb
		    flag = 0;
		    relres = normr0;
		    break;
	    end
        % Step
        s = M \ q;
        b_top = r' * s;
        b_bot = q' * s; % breakdown if b_bot is zero
        lambda = b_top / b_bot;
    
        x = x + lambda*p;
        r = r - lambda*q;

        feasible = true;
        for i = 1:length(x)
            if(mod(mindices(i),3)==1 && x(i) < 1e-6)
                feasible = false;
            end

            if(mod(mindices(i),3)==2 && norm([x(i) x(i+1)]) > mu*x(i-1)+ 1e-6)
                feasible = false;
            end
        end
    
        if(~feasible)
            flag = 1;
            break;
        end
    
        t = A * s;
    
        alpha0 = (t' * s) / b_bot;
        alpha1 = (t' * s_old) / b_bot_old;
    
        p_temp = p;
        q_temp = q;
    
        p = s - alpha0 * p - alpha1 * p_old;
        q = t - alpha0 * q - alpha1 * q_old;
    
        s_old = s;
        p_old = p_temp;
        q_old = q_temp;
        b_bot_old = b_bot;
    end
end

%% Gauss-Seidel 
function [x,flag,relres,iter,resvec] = gs(A,b,tol,maxit,M,~,x0,~,uReduced,~,mindices, mu)
    flag = 1;
    relres = 0;
    if nargout >= 5
	    resvec = [];
    end
    
    x = x0;
    r = b - A * x;
    normb = norm(b);
    for iter = 1 : maxit
	    % Convergence check
	    normr0 = norm(r);
	    if nargout >= 5
		    resvec(end+1) = normr0; %#ok<AGROW>
	    end
	    if normr0 < tol*normb
		    flag = 0;
		    relres = normr0;
		    break;
        end

        feasible = true;
        for i = 1 : length(x)
            ri = b(i) - A(i,:)*x;
            x(i) = x(i) + ri / A(i,i);
            if(mod(mindices(i),3)==1 && x(i) < 0)
                feasible = false;
            end
            if(mod(mindices(i),3)==2 && norm([x(i) x(i+1)]) > mu*x(i-1))
                feasible = false;
            end
        end

        if(~feasible)
            flag = 1;
            break;
        end
        r = b - A * x;
    end
end

function [xc, xd, tIndex]= computeCauchyPointCone(G,c,l,u,x,mu)
    g = G * x + c;
    n = size(g, 1);
    tList = Inf(n, 1);
    xt = zeros(n, 1);
    xc = xt;
    for i = 1:3:n
        t = rayConeIntersection(x(i:i+2),-g(i:i+2),mu);
        tList(i) = t;
        if((abs(g(i)) < 1e-6 && abs(x(i))<1e-6))
            tList(i+1) = Inf;
        elseif(g(i)>0)
            tList(i+1) = x(i) / g(i);
        else
            tList(i+1) = Inf;
        end
        tList(i+2) = Inf;
    end

    tUniqueList = unique(tList,'sorted');
    if(tUniqueList(end) ~= Inf)
        tUniqueList(end + 1) = Inf;
    end
    t = 0;
    for i = 1: size(tUniqueList,1)
        %{
        tIndex = t < tList;
        temp = x - t * g;
        xt(tIndex) = temp(tIndex);
        temp = x - tList .* g;
        xt(~tIndex) = temp(~tIndex);
        p = zeros(n, 1);
        p(tIndex) = -g(tIndex);
        %}
        [xt,p] = compute_xt_p(tList, t, x, g);
        fPrime = c' * p + xt' * G * p;
        fPrimePrime = p' * G * p;
        deltaTStar = - fPrime / fPrimePrime;
        if(fPrime >0)
            break;
        elseif(deltaTStar >=0 && deltaTStar < tUniqueList(i) - t)
            t = t + deltaTStar;
            break;
        end
        t = tUniqueList(i);
    end
    %{
    tIndex = t < tList;
    temp = x - t * g;
    xc(tIndex) = temp(tIndex);
    temp = x - tList .* g;
    xc(~tIndex) = temp(~tIndex);
    %}
    [xc,~] = compute_xt_p(tList, t, x, g);

    tIndex = false(n,1);
    xd = zeros(n,1);
    for i = 1:3:n
        if(t<tList(i))
            tIndex(i:i+2) = true;
        elseif(t<tList(i+1))
            tIndex(i) = true;
        else
            tIndex(i:i+2) = false;
        end
        if(norm(xc(i:i+2))<1e-9)
            xd(i:i+2) = [1 0 0]';
        else
            xd(i:i+2) = xc(i:i+2) / norm(xc(i:i+2));
        end
    end
end

function [xt,p] = compute_xt_p(tList,t,x,g)
    n = length(x);
    xt = zeros(n,1);
    p = zeros(n,1);
    for i = 1:3:n
        if(t < tList(i))
            xt(i:i+2) = x(i:i+2) - t*g(i:i+2);
            p(i:i+2) = -g(i:i+2);
        elseif(t < tList(i+1))
            gp = g(i:i+2);
            if(norm(x(i:i+2))<1e-9)
                xd = [1 0 0]';
            else
                xd = x(i:i+2) / norm(x(i:i+2));
            end
            gp = (xd'*gp)*xd;
            xt(i:i+2) = x(i:i+2) - tList(i)*g(i:i+2) - (t - tList(i))*gp;
            p(i:i+2) = -gp;
        else
            gp = g(i:i+2);
            if(norm(x(i:i+2))<1e-9)
                xd = [1 0 0]';
            else
                xd = x(i:i+2) / norm(x(i:i+2));
            end
            gp = (xd'*gp)*xd;
            xt(i:i+2) = x(i:i+2) - tList(i)*g(i:i+2) - (tList(i+1) - tList(i))*gp;
            p(i:i+2) = zeros(3,1);
        end
    end

end

%% Ray-Cone Intersection
function t = rayConeIntersection(x, g, mu)
    gnorm = norm(g);
    g = g / gnorm;

    A = g(2)^2 + g(3)^2 - g(1)^2 * mu^2;
    B = 2 * (x(2)*g(2) + x(3)*g(3) - x(1)*g(1) * mu^2);
    C = x(2)^2 + x(3)^2 - x(1)^2 * mu^2;

    if(abs(C)<1e-6 && gnorm < 1e-6)
        t = Inf;
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
        t = t / gnorm;
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
        t = t / gnorm;
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
        t = t / gnorm;
        return;
    end

    if(t2 ==0)
        if(x(1) + t1*g(1) >0)
            t = t2;
        else
            t = Inf;
        end
        t = t / gnorm;
        return;
    end
    
    
    if(t1>0 && t2 >0)
        t = min(t1,t2);
    elseif(t1<0 && t2<0)
        t = Inf;
    else
        t = max(t1,t2);
    end
    t = t / gnorm;
end
