mu = 0.5;
load("test_data.mat")
options.ProjectionMethod = 'direct';
options.MaxIterations = 1000;
options.CGMaxIterations=200;
options.Tolerance = 1e-9;

l = -inf(length(b),1);
for i = 1:3:length(b)
    l(i) = 0;
end
u = inf(length(b),1);
x = zeros(length(b),1);
[x, f, exitflag, output, lambda]= cone_gpqp(A,-b,l,u,x,options,mu);