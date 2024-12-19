%% Numerical Optimization (Jorge Nocedal and Stephen J. Wright), Springer, 2006.
%% Gradient Projection Method for QP (Algorithm 16.5)
function [x, f, exitflag, output, lambda]= gpqp(G,c,l,u,x,opts)
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

fPrev = 0.5 * x' * G * x + x' * c;
gPrev = G * x + c;
CGiterVec = [];
neqIndexPrev = zeros(n);
for i = 1:iterNum
    [xc,neqIndex] = computeCauchyPoint(G,c,l,u,x);
    lNeq = l(neqIndex);
    uNeq = u(neqIndex);
    %xY = A * xc;
    xY = xc(~neqIndex);
    %xZ = Z' * xc;
    xZ = xc(neqIndex);
    %cZ = Z' * (G * (Y * xY)) + Z' * c;
    if(isempty(xY))
        cZ = c(neqIndex);
    else
        cZ = G(neqIndex, ~neqIndex) * xY + c(neqIndex);
    end

    %H = eye(n);
    %H = diag(abs(diag(G)));
    %W_zz = Z' * H * Z;
    H = abs(diag(G));
    W_zz = H(neqIndex);

    %xZ = pcg_rs(Z' * G * Z, -cZ,eps,CGiterNum,W_zz,[],xZ,lNeq,uNeq,projectionOpts);
    [xZ,~,~,CGiter,~] = pcg_rs(G(neqIndex, neqIndex), -cZ,eps,CGiterNum,W_zz,[],xZ,lNeq,uNeq,projectionOpts);
    CGiterVec = [CGiterVec, CGiter];

    x(~neqIndex) = xc(~neqIndex);
    x(neqIndex) = xZ;

    % if satisfies the KKT conditions
    f = 0.5 * x' * G * x + x' * c;
    g = G * x + c;

	if norm(g) < eps && sum(neqIndex) == n
		break;
	end
    if norm(g - gPrev) < eps
        break;
    end
    if norm(neqIndex - neqIndexPrev) < eps
        %break;
    end
    fPrev = f;
    gPrev = g;
    neqIndexPrev = neqIndex;
end

lambdaLower = G * x + c;
lambdaUpper = lambdaLower;
lambdaLower(lambdaLower < 0) = 0;
lambdaUpper(lambdaUpper > 0) = 0;
lambda.lower = lambdaLower;
lambda.upper = -lambdaUpper;

if(i ~= iterNum)
    exitflag = 1;
else
    exitflag = 0;
end

output.iterations = i;
output.projectionMethod = projectionOpts;
output.cgiterations = CGiterVec;
end

%% Preconditioned CG for Reduced Systems (Algorithm 16.1 )
function [xZ, flag, relres, j, resvec] = pcg_rs(A,b,tol,maxit,M1,~,x0,lReduced,uReduced,projectionOpts)
%M1 is assumed to be a diagnal matrix, so it is represented by a vector
%inside this function.
n = size(b,1);
if isempty(tol)
    tol = 1e-9;
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
    if(norm(rZ' * M1 * rZ) < tol)
        flag = 0;
        break;
    end
    alpha = (rZ' * gZ) / (dZ' * A * dZ);
    xZ_last = xZ;
    xZ = xZ + alpha * dZ;

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

%% Computes Cauchy Point (Algorithm 16.5 line 5)
function [xc, tIndex]= computeCauchyPoint(G,c,l,u,x)
g = G * x + c;
n = size(g, 1);
tList = Inf(n, 1);
xt = zeros(n, 1);
xc = xt;
lBoundIndex = (g>0) & (l > -Inf );
uBoundIndex = (g<0) & (u < Inf );
temp = (x - u) ./ g;
tList(uBoundIndex) = temp(uBoundIndex);
temp = (x - l) ./ g;
tList(lBoundIndex) = temp(lBoundIndex);
tUniqueList = unique(tList,'sorted');
if(tUniqueList(end) ~= Inf)
    tUniqueList(end + 1) = Inf;
end
t = 0;
for i = 1: size(tUniqueList,1)
    tIndex = t <= tList;
    temp = x - t * g;
    xt(tIndex) = temp(tIndex);
    temp = x - tList .* g;
    xt(~tIndex) = temp(~tIndex);
    p = zeros(n, 1);
    p(tIndex) = -g(tIndex);
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
tIndex = t < tList;
temp = x - t * g;
xc(tIndex) = temp(tIndex);
temp = x - tList .* g;
xc(~tIndex) = temp(~tIndex);
end