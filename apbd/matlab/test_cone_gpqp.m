%{
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
%}

x = [2.44999999992264
-1.01070349431031e-10
-7.94194419907315e-11];
g = -[4.26325641456060e-13
-8.09368754905583e-10
-1.06788792015864e-09]';
mu = 0.5;
t = rayConeIntersection(x, g, mu);

%% Ray-Cone Intersection
function t = rayConeIntersection(x, g, mu)
    % rayConeIntersection calculates the intersection of a ray with a cone.
    %
    % Inputs:
    %   x         - 3D start point of the ray [x1; x2; x3]
    %   g         - 3D direction vector of the ray [g1; g2; g3]
    %   mu        - Fricitonal coefficient
    %
    % Output:
    %   t         - Intersection parameter(s). Empty if no intersection.


    % Coefficients of the quadratic equation (At^2 + Bt + C = 0)
    gnorm = norm(g);
    g = g / gnorm;

    A = g(2)^2 + g(3)^2 - g(1)^2 * mu^2;
    B = 2 * (x(2)*g(2) + x(3)*g(3) - x(1)*g(1) * mu^2);
    C = x(2)^2 + x(3)^2 - x(1)^2 * mu^2;

    if(abs(C)<1e-9)
        C = 0;
    end

    if(abs(C)<1e-9 && gnorm < 1e-9)
        t = 0;
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