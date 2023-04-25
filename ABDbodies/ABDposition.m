clear;
rng(2);

h = 0.01;
tEnd = 3;
drawHz = 60;

%E = eye(4);
E = se3.randE();
%E(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
E(1:3,4) = [0 0 4]';

% Note: whenever we use reshape(), don't forget that we want q to be
% row-wise, so we need the transpose.
q = [E(1:3,4); reshape(E(1:3,1:3)',9,1)];
phi = [10 1 0 0 0 0]'; % rigid velocity in local space
Edot = E*se3.brac(phi);
qdot = zeros(12,1);
qdot(1:3) = Edot(1:3,4);
qdot(4:12) = reshape(Edot(1:3,1:3)',9,1);

grav = [0 0 -9.8]';
sides = [2 1 0.5];
density = 1;
volume = prod(sides);
mass = volume*density;

% I is the 6x1 diagonal rigid inertia, assuming that the frame origin is at the
% center of mass and the axes are oriented along the principal axes. We store
% the rotation on top of translations, so that I(1:3) is the rotational inertia
% and I(4:6) is the translational inertia.
I = se3.inertiaCuboid(sides,density);
% This matrix transforms from rigid to affine.
R2A = [
	-1  1  1
	 1 -1  1
	 1  1 -1
	]/2;
% Transform translational inertia
M(1:3) = I(4:6);
% Transform rotational inertia, repeating 3 times
M(4:6) = R2A*I(1:3);
M(7:9) = M(4:6);
M(10:12) = M(4:6);
M = M';

% Sim loop
for t = 0 : h : tEnd
	q0 = q;
	% Unconstrained
	f = zeros(12,1);
	f(1:3) = M(1:3).*grav;
	q = q + h*qdot + h^2*(f./M);
	% Volume
	[C,dC] = ABDvolume(q);
	lambda = -C/(dC*(dC'./M));
	q = q + lambda*(dC'./M);
	% Ortho
	[C,dC] = ABDortho(q);
	for i = 1 : length(C)
		lambda = -C(i)/(dC(i,:)*(dC(i,:)'./M));
		q = q + lambda*(dC(i,:)'./M);
	end
	% Collisions
	collisions = floorCollision(q,sides);
	for i = 1 : length(collisions)
		c = collisions(i);
		C = c.nw'*c.xw;
		dC = c.nw'*ABDjacobian(c.xl);
		lambda = -C/(dC*(dC'./M));
		q = q + lambda*(dC'./M);
	end
	qdot = (q - q0)/h;
	% Draw
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		draw(t,q,sides);
	end
end

%%
function J = ABDjacobian(xl)
J = zeros(3,12);
J(1:3,1:3) = eye(3);
J(1,4:6) = xl';
J(2,7:9) = xl';
J(3,10:12) = xl';
end

%%
function [C,dC] = ABDvolume(q)
% clear
% syms A [3,3] real
% C = det(A) - 1;
% dC(4:12) = jacobian(C,reshape(A',9,1));
A = reshape(q(4:12),3,3)';
C = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) - 1;
dC = [0, 0, 0, A(2,2)*A(3,3) - A(2,3)*A(3,2), A(2,3)*A(3,1) - A(2,1)*A(3,3), A(2,1)*A(3,2) - A(2,2)*A(3,1), A(1,3)*A(3,2) - A(1,2)*A(3,3), A(1,1)*A(3,3) - A(1,3)*A(3,1), A(1,2)*A(3,1) - A(1,1)*A(3,2), A(1,2)*A(2,3) - A(1,3)*A(2,2), A(1,3)*A(2,1) - A(1,1)*A(2,3), A(1,1)*A(2,2) - A(1,2)*A(2,1)];
end

%%
function [C,dC] = ABDortho(q)
% clear
% syms A [3,3] real
% C(1,1) = A(:,1)'*A(:,1) - 1;
% C(2,1) = A(:,2)'*A(:,2) - 1;
% C(3,1) = A(:,3)'*A(:,3) - 1;
% C(4,1) = A(:,1)'*A(:,2);
% C(5,1) = A(:,2)'*A(:,3);
% C(6,1) = A(:,3)'*A(:,1);
% dC(:,4:12) = jacobian(C,reshape(A',9,1));
A = reshape(q(4:12),3,3)';
C(1,1) = A(1,1)^2 + A(2,1)^2 + A(3,1)^2 - 1;
C(2,1) = A(1,2)^2 + A(2,2)^2 + A(3,2)^2 - 1;
C(3,1) = A(1,3)^2 + A(2,3)^2 + A(3,3)^2 - 1;
C(4,1) = A(1,1)*A(1,2) + A(2,1)*A(2,2) + A(3,1)*A(3,2);
C(5,1) = A(1,2)*A(1,3) + A(2,2)*A(2,3) + A(3,2)*A(3,3);
C(6,1) = A(1,1)*A(1,3) + A(2,1)*A(2,3) + A(3,1)*A(3,3); 
dC = [
    0, 0, 0, 2*A(1,1),        0,        0, 2*A(2,1),        0,        0, 2*A(3,1),        0,        0
    0, 0, 0,        0, 2*A(1,2),        0,        0, 2*A(2,2),        0,        0, 2*A(3,2),        0
    0, 0, 0,        0,        0, 2*A(1,3),        0,        0, 2*A(2,3),        0,        0, 2*A(3,3)
    0, 0, 0,   A(1,2),   A(1,1),        0,   A(2,2),   A(2,1),        0,   A(3,2),   A(3,1),        0
    0, 0, 0,        0,   A(1,3),   A(1,2),        0,   A(2,3),   A(2,2),        0,   A(3,3),   A(3,2)
    0, 0, 0,   A(1,3),        0,   A(1,1),   A(2,3),        0,   A(2,1),   A(3,3),        0,   A(3,1)
];
end

%%
function [C,dC] = ABDortho1(q)
% clear
% syms A [3,3] real
% C = reshape(A*A' - eye(3),9,1);
% dC(:,4:12) = jacobian(C,reshape(A',9,1));
A = reshape(q(4:12),3,3)';
C = [
	A(1,1)^2 + A(1,2)^2 + A(1,3)^2 - 1
	A(1,1)*A(2,1) + A(1,2)*A(2,2) + A(1,3)*A(2,3)
	A(1,1)*A(3,1) + A(1,2)*A(3,2) + A(1,3)*A(3,3)
	A(1,1)*A(2,1) + A(1,2)*A(2,2) + A(1,3)*A(2,3)
	A(2,1)^2 + A(2,2)^2 + A(2,3)^2 - 1
	A(2,1)*A(3,1) + A(2,2)*A(3,2) + A(2,3)*A(3,3)
	A(1,1)*A(3,1) + A(1,2)*A(3,2) + A(1,3)*A(3,3)
	A(2,1)*A(3,1) + A(2,2)*A(3,2) + A(2,3)*A(3,3)
	A(3,1)^2 + A(3,2)^2 + A(3,3)^2 - 1
	];
dC = [
	0, 0, 0, 2*A(1,1), 2*A(1,2), 2*A(1,3),        0,        0,        0,        0,        0,        0
	0, 0, 0,   A(2,1),   A(2,2),   A(2,3),   A(1,1),   A(1,2),   A(1,3),        0,        0,        0
	0, 0, 0,   A(3,1),   A(3,2),   A(3,3),        0,        0,        0,   A(1,1),   A(1,2),   A(1,3)
	0, 0, 0,   A(2,1),   A(2,2),   A(2,3),   A(1,1),   A(1,2),   A(1,3),        0,        0,        0
	0, 0, 0,        0,        0,        0, 2*A(2,1), 2*A(2,2), 2*A(2,3),        0,        0,        0
	0, 0, 0,        0,        0,        0,   A(3,1),   A(3,2),   A(3,3),   A(2,1),   A(2,2),   A(2,3)
	0, 0, 0,   A(3,1),   A(3,2),   A(3,3),        0,        0,        0,   A(1,1),   A(1,2),   A(1,3)
	0, 0, 0,        0,        0,        0,   A(3,1),   A(3,2),   A(3,3),   A(2,1),   A(2,2),   A(2,3)
	0, 0, 0,        0,        0,        0,        0,        0,        0, 2*A(3,1), 2*A(3,2), 2*A(3,3)
	];
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
