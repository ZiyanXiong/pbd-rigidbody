load("layer3matrix.mat");
A=Al;
b = bl;
n = length(b);
x = zeros(n,1);
dx = zeros(n,1);
upiterMax = 150;
mu = 0.5;

inds = 1:n;
nind = inds(1:3:end);
tind = inds(2:3:end);
for iter = 1:upiterMax
    r0 = b-A*dx;

    for i = nind
        ri = b(i) - A(i,:) * dx;
        dx(i) = dx(i) + ri ./ A(i,i); 
        
        if(mod(i,3) == 1)
            if(x(i) + dx(i) < 0)
                dx(i) = -x(i);
            end
        end
    end

    for i = tind
        ri = b(i:i+1) - A(i:i+1,:) * dx;
        dx(i:i+1) = dx(i:i+1) + A(i:i+1,i:i+1)\ri; 
        
        if (norm([x(i) + dx(i) x(i+1) + dx(i+1)]) > mu * (x(i-1)+dx(i-1)))
            scale =  mu * (x(i-1)+dx(i-1)) / norm([x(i) + dx(i) x(i+1) + dx(i+1)]);
            dx(i) = scale * (x(i) + dx(i)) - x(i);
            dx(i+1) = scale * (x(i+1) + dx(i+1)) - x(i+1);
        end
    end

    r = b-A*dx;
    if(norm(r(1:3:n)-r0(1:3:n)) < 1e-4)
        break;
    end
end
disp(iter)