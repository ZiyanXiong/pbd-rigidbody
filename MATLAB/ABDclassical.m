clear;
rng(0);

h = 0.01;
tEnd = 3;

E = se3.randE();
R = E(1:3,1:3);

q = [[0 0 4]'; reshape(R,9,1)];
qdot = zeros(12,1);

grav = [0 0 -9.8]';
sides = [2 1 0.5];
density = 1;
volume = prod(sides);
mass = volume*density;

stiffness = 1e3;
floorStiffness = 1e3;
damping = 1e-1;
D = damping*eye(12);

% Sample the cuboid to build the mass matrix
% To double check: if the body is centered, then we get a diagonal matrix
nx = 10;
ny = 10;
nz = 10;
n = nx*ny*nz;
m = mass/n;
M = zeros(12,12);
for sx = linspace(-0.5*sides(1),0.5*sides(1),nx)
	for sy = linspace(-0.5*sides(2),0.5*sides(2),ny)
		for sz = linspace(-0.5*sides(3),0.5*sides(3),nz)
			xl = [sx sy sz]';
			J = ABDjacobian(xl);
			M = M + m*(J'*J);
		end
	end
end

% Sim loop
for t = 0 : h : tEnd
	f = zeros(12,1);
	f(1:3) = f(1:3) + mass*grav;
	collisions = floorCollision(q,sides);
	for i = 1 : length(collisions)
		c = collisions(i);
		Jcol = ABDjacobian(c.xl);
		fcol = -floorStiffness*c.d*c.nw;
		f = f + Jcol'*fcol;
	end
	[~,gorth,Korth] = ABDorth(q(4:end));
	f(4:end) = f(4:end) - stiffness*gorth;
	K = zeros(12,12);
	K(4:end,4:end) = -stiffness*Korth;
	LHS = M + h*D - h^2*K;
	rhs = M*qdot + h*f;
	qdot = LHS\rhs;
	q = q + h*qdot;
	draw(t,q,sides);
end

%%
function J = ABDjacobian(xl)
J = zeros(3,12);
J(1:3,1:3) = eye(3);
J(1,4:6) = xl';
J(2,7:9) = xl';
J(3,10:12) = xl';
end

function [V,g,H] = ABDorth(A)
%A = reshape(A,3,3);
%V = norm(A*A' - eye(3),'fro')^2;
% Note that these three vectors are column vectors, but they represent the
% three rows of A.
a1 = A(1:3);
a2 = A(4:6);
a3 = A(7:9);
% Inner products
a11 = a1'*a1;
a12 = a1'*a2;
a13 = a1'*a3;
a21 = a12;
a22 = a2'*a2;
a23 = a2'*a3;
a31 = a13;
a32 = a23;
a33 = a3'*a3;
% Potential
V = (a11 - 1)^2 + (a22 - 1)^2 + (a33 - 1)^2 ...
	                + a12^2 + a13^2 ...
	+ a21^2                 + a23^2 ...
	+ a31^2 + a32^2                ;
if nargout >= 2
	% Gradient
	g = zeros(9,1);
	g(1:3) = 4*((a11 - 1)*a1 + a12*a2 + a13*a3);
	g(4:6) = 4*((a22 - 1)*a2 + a23*a3 + a21*a1);
	g(7:9) = 4*((a33 - 1)*a3 + a31*a1 + a32*a2);
	if nargout >= 3
		% Hessian
		H = zeros(9,9);
		I = eye(3);
		% Outer products
		A11 = a1*a1';
		A12 = a1*a2';
		A13 = a1*a3';
		A21 = A12';
		A22 = a2*a2';
		A23 = a2*a3';
		A31 = A13';
		A32 = A23';
		A33 = a3*a3';
		H(1:3,1:3) = 4*((2*A11 + (a11 - 1)*I) + A22 + A33);
		H(1:3,4:6) = 4*(A21 + a12*I);
		H(1:3,7:9) = 4*(A31 + a13*I);
		H(4:6,1:3) = H(1:3,4:6)';
		H(4:6,4:6) = 4*((2*A22 + (a22 - 1)*I) + A33 + A11);
		H(4:6,7:9) = 4*(A32 + a23*I);
		H(7:9,1:3) = H(1:3,7:9)';
		H(7:9,4:6) = H(4:6,7:9)';
		H(7:9,7:9) = 4*((2*A33 + (a33 - 1)*I) + A11 + A22);
	end
end
end

%%
function collisions = floorCollision(q,sides)
E = eye(4);
E(1:3,4) = q(1:3);
E(1:3,1:3) = reshape(q(4:12),3,3)'; % row-by-row
sides = sides/2;
xl = ones(4,8);
xl(1:3,1) = [-sides(1), -sides(2), -sides(3)]';
xl(1:3,2) = [ sides(1), -sides(2), -sides(3)]';
xl(1:3,3) = [ sides(1),  sides(2), -sides(3)]';
xl(1:3,4) = [-sides(1),  sides(2), -sides(3)]';
xl(1:3,5) = [-sides(1), -sides(2),  sides(3)]';
xl(1:3,6) = [ sides(1), -sides(2),  sides(3)]';
xl(1:3,7) = [ sides(1),  sides(2),  sides(3)]';
xl(1:3,8) = [-sides(1),  sides(2),  sides(3)]';
xw = E*xl;
collisions = [];
for i = 1 : 8
	if xw(3,i) < 0
		collisions(end+1).xl = xl(1:3,i); %#ok<AGROW>
		collisions(end).xw = xw(1:3,i);
		collisions(end).nw = [0 0 1]';
		collisions(end).d = xw(3,i); % penetration depth
	end
end
end

%%
function draw(t,q,sides)
if t == 0
	clf;
	hold on;
	axis equal;
	axis(3*[-1 1 -1 1 0 2]);
	grid on;
	view(3);
	xlabel('X');
	ylabel('Y');
	zlabel('Z');
	ax = gca;
	ax.Clipping = 'off';
end

cla;
E = eye(4);
E(1:3,4) = q(1:3);
E(1:3,1:3) = reshape(q(4:12),3,3)'; % row-by-row
se3.drawAxis(E);
se3.drawCuboid(E,sides);
title(sprintf('%f',t));
drawnow
end
