%% Numerical Optimization (Jorge Nocedal and Stephen J. Wright), Springer, 2006.
%% Gradient Projection Method for QP (Algorithm 16.5)
function [x, f, exitflag, output]= matrixfree_gpqp(G,c,l,u,x,opts, mu)
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

rs = zeros(iterNum,1);
g = zeros(n,1);
x = zeros(n,1);
nactiveIndex = true(n,1);
blockList = {};
for i = 1:3:n
    blockList{end+1} = GPQPBlock(G,c,x,i);
end

fPrev = 0;
gPrev = g;
for i = 1: length(blockList)
    blockList{i}.blockList = blockList;
    blockList{i}.compute_Ax_g();
    fPrev = fPrev + blockList{i}.x' * ( 0.5 * blockList{i}.Ax +  blockList{i}.b);
    gPrev(blockList{i}.ind:blockList{i}.ind+2) = blockList{i}.g;
end
CGiterVec = [];

for iter = 1:iterNum
    if(iter == 3)
        disp(iter);
    end
    computeCauchyPoint(blockList, n);

    for i = 1: length(blockList)
        if (norm([blockList{i}.x(2:3)]) > mu * blockList{i}.x(1))
            blockList{i}.activeInd(2:3) = true;
        end
        blockList{i}.l(2:3) = -mu*blockList{i}.x(1);
        blockList{i}.u(2:3) = mu*blockList{i}.x(1);

        blockList{i}.compute_Ax_g();
        g(blockList{i}.ind:blockList{i}.ind+2) = blockList{i}.g;
        x(blockList{i}.ind:blockList{i}.ind+2) = blockList{i}.x;
        nactiveIndex(blockList{i}.ind:blockList{i}.ind+2) = ~blockList{i}.activeInd;
    end

    [~,CGiter,~] = pcg_rs(blockList, [], []);

    CGiterVec = [CGiterVec, CGiter];

    for j = 1: length(blockList)
        blockList{j}.project(mu);
    end

    % if satisfies the KKT conditions
    f = 0;
    for i = 1: length(blockList)
        blockList{i}.compute_Ax_g();
        f = f + blockList{i}.x' * ( 0.5 * blockList{i}.Ax +  blockList{i}.b);
        g(blockList{i}.ind:blockList{i}.ind+2) = blockList{i}.g;
        x(blockList{i}.ind:blockList{i}.ind+2) = blockList{i}.x;
    end

	if norm(g) < eps
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
    rs(iter) = norm(g);
end

if(i ~= iterNum)
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
function [flag, CGiters, resvec] = pcg_rs(blockList,tol,maxit)
    %M1 is assumed to be a diagnal matrix, so it is represented by a vector
    %inside this function.

    if isempty(tol)
        tol = 1e-9;
    end
    if isempty(maxit)
        maxit = 100;
    end
    for i = 1: length(blockList)
        blockList{i}.Minv_cg = ones(3,1); %Preconditioner
        blockList{i}.compute_b_cg();
        blockList{i}.init_cg() ;
    end
    
    resvec = [];
    for CGiters = 1: maxit
        r = 0;
        for i = 1: length(blockList)
            M = 1 ./ blockList{i}.Minv_cg;
            r = r + blockList{i}.r_cg(~blockList{i}.activeInd)' * diag(M(~blockList{i}.activeInd)) * blockList{i}.r_cg(~blockList{i}.activeInd);
        end
    
        if(r < tol)
            break;
        end
        numerator = 0;
        denominator = 0;
        for i = 1: length(blockList)
            blockList{i}.compute_Ad_cg();
            numerator = numerator + blockList{i}.r_cg(~blockList{i}.activeInd)' * blockList{i}.g_cg(~blockList{i}.activeInd);
            denominator = denominator + blockList{i}.d_cg(~blockList{i}.activeInd)' * blockList{i}.Ax(~blockList{i}.activeInd);
        end
        alpha = numerator / denominator;
    
        feasible = true;
        for i = 1: length(blockList)
            feasible = feasible & blockList{i}.update_cg(alpha);
        end
    
        if(~feasible)
            break;
        end
        numerator = 0;
        denominator = 0;
        for i = 1: length(blockList)
            denominator = denominator + blockList{i}.r_cg(~blockList{i}.activeInd)' *  blockList{i}.g_cg(~blockList{i}.activeInd);
            blockList{i}.r_cg = blockList{i}.r_cg + alpha * blockList{i}.Ax;
            blockList{i}.g_cg = blockList{i}.Minv_cg .* blockList{i}.r_cg;
            numerator = numerator + blockList{i}.r_cg(~blockList{i}.activeInd)' *  blockList{i}.g_cg(~blockList{i}.activeInd);
        end
        beta = numerator / denominator;
        for i = 1: length(blockList)
            blockList{i}.d_cg = -blockList{i}.g_cg + beta * blockList{i}.d_cg;
        end
    end
    if(CGiters == maxit)
        flag = 2;
    else
        flag = 1;
    end
end

%% Computes Cauchy Point (Algorithm 16.5 line 5)
function computeCauchyPoint(blockList,n)
    tList = Inf(n, 1);
    for i = 1: length(blockList)
        blockList{i}.compute_Ax_g();
        tList(blockList{i}.ind:blockList{i}.ind+2) = blockList{i}.compute_tbar();
    end
    tUniqueList = unique(tList,'sorted');
    if(tUniqueList(end) ~= Inf)
        tUniqueList(end + 1) = Inf;
    end
    t = 0;
    for i = 1: size(tUniqueList,1)
        for j = 1:length(blockList)
            blockList{j}.compute_p(t);
            blockList{j}.compute_xc(t);
        end
        fPrime = 0;
        fPrimePrime = 0;
        for j = 1:length(blockList)
            blockList{j}.compute_Gp();
            fPrime = fPrime + blockList{j}.b' * blockList{j}.p +  blockList{j}.xc' * blockList{j}.Gp;
            fPrimePrime = fPrimePrime + blockList{j}.p' * blockList{j}.Gp;
        end
        deltaTStar = - fPrime / fPrimePrime;
        if(fPrime >0)
            break;
        elseif(deltaTStar >=0 && deltaTStar < tUniqueList(i) - t)
            t = t + deltaTStar;
            break;
        end
        t = tUniqueList(i);
    end
    for i = 1: length(blockList)
        blockList{i}.compute_xc(t);
        blockList{i}.x = blockList{i}.xc;
    end
end