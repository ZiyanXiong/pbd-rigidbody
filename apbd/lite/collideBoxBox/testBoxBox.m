% Test driver for collideBoxBox

if false
	s1 = [1 1 1];
	s2 = [2 2 2];
	E1 = eye(4);
	E2 = [eye(3),[1 1 -1.45]'; 0 0 0 1];
else
	% Good cases: 1 5 9 19
	rng(5);
	s1 = [3 2 2];
	s2 = [3 2 1];
	E1 = randE();
	E2 = randE();
end

% flip 1 and 2
if false
	s0 = s1;
	s1 = s2;
	s2 = s0;
	E0 = E1;
	E1 = E2;
	E2 = E0;
end

clf;
hold on;
axis equal;
view(3);
grid on;
[F1,V1] = patchCuboid(E1,s1);
[F2,V2] = patchCuboid(E2,s2);
patch('Faces',F1,'Vertices',V1,'FaceColor',[1.0 0.5 0.5]);
patch('Faces',F2,'Vertices',V2,'FaceColor',[0.5 1.0 0.5]);
xlabel('X');
ylabel('Y');
zlabel('Z');
alpha(0.5);

thresh = 1e-6;
R1 = E1(1:3,1:3);
R2 = E2(1:3,1:3);
collisions = odeBoxBox_mex(E1,s1,E2,s2);
nw = collisions.nor; % The normal is outward from body 1 (red)
n1 =  R1'*nw;
n2 = -R2'*nw; % negate since nw is defined wrt body 1
bmin1 = -0.5*s1';
bmax1 =  0.5*s1';
bmin2 = -0.5*s2';
bmax2 =  0.5*s2';
for i = 1 : collisions.count
	xw = collisions.pos(:,i);
	d = collisions.depth(i);
	xnw = [xw,xw+nw];
	% Plot position and normal in world coords
	plot3(xw(1),xw(2),xw(3),'kx');
	plot3(xnw(1,:),xnw(2,:),xnw(3,:),'r-');
	% Compute local point on body 1 with ray casting
	x1 = E1\[xw;1];
	x1 = x1(1:3);
	x1 = (1 - thresh)*x1; % make the point go slightly inside the box
	[hit,t] = smits_mul(x1,-n1,bmin1,bmax1); % negate ray since it starts inside the box
	assert(hit); % There should be a hit since the ray started inside the box
	x1 = x1 - t*n1; % negate since smits_mul returns negative t for rays starting inside the box
	x1w = E1*[x1;1];
	plot3(x1w(1),x1w(2),x1w(3),'ro');
	% Compute local point on body 2 with ray casting
	x2 = E2\[xw;1];
	x2 = x2(1:3);
	x2 = (1 - thresh)*x2; % make the point go slightly inside the box
	[hit,t] = smits_mul(x2,-n2,bmin2,bmax2); % negate ray since it starts inside the box
	assert(hit); % There should be a hit since the ray started inside the box
	x2 = x2 - t*n2; % negate since smits_mul returns negative t for rays starting inside the box
	x2w = E2*[x2;1];
	plot3(x2w(1),x2w(2),x2w(3),'go');
end

%%

clear
s1 = [3 2 2];
E1 = eye(4);
clf
hold on
axis equal
view(3)
grid on
[F1,V1] = patchCuboid(E1,s1);
patch('Faces',F1,'Vertices',V1,'FaceColor',[1.0 0.5 0.5]);
alpha(0.5);
bmin1 = -s1'/2;
bmax1 =  s1'/2;
rx = [0 0 0]';
for k = 1 : 20
	rv = randn(3,1);
	rv = rv/norm(rv);
	[hit,t] = smits_mul(rx, rv, bmin1, bmax1);
	t = -t; % if the origin is inside the box, the result is negative
	if hit
		x = rx + t*rv;
		rxv = [rx,rx+norm(bmax1)*rv];
		plot3(rxv(1,:),rxv(2,:),rxv(3,:),'k-');
		plot3(x(1),x(2),x(3),'rx');
	end
end

%%

clear;
data = load('../data.mat');
rx = 0.9999*data.rx;
rv = data.rv;
bmin = data.bmin;
bmax = data.bmax;
s1 = -2*bmin';
E1 = eye(4);
clf
hold on
axis equal
view(3)
grid on
[F1,V1] = patchCuboid(E1,s1);
patch('Faces',F1,'Vertices',V1,'FaceColor',[1.0 0.5 0.5]);
alpha(0.5);

plot3(rx(1),rx(2),rx(3),'rx');
rxv = [rx,rx+norm(bmax)*rv];
plot3(rxv(1,:),rxv(2,:),rxv(3,:),'k-');
[hit,t] = smits_mul(rx, rv, bmin, bmax);

t = -t; % if the origin is inside the box, the result is negative
if hit
	x = rx + t*rv;
	rxv = [rx,rx+norm(bmax)*rv];
	plot3(rxv(1,:),rxv(2,:),rxv(3,:),'k-');
	plot3(x(1),x(2),x(3),'rx');
end

%%
function E = randE()
% Creates a random transformation matrix
[Q,R] = qr(randn(3)); %#ok<ASGLU>
if det(Q) < 0
	% Negate the Z-axis
	Q(:,3) = -Q(:,3);
end
E = [Q, randn(3,1); 0 0 0 1];
end

%%
function [F,V] = patchCuboid(E,whd)
% Gets the faces and vertices for a cuboid at E with (width, height, depth)
whd = whd/2;
verts = ones(4,8);
verts(1:3,1) = [-whd(1), -whd(2), -whd(3)]';
verts(1:3,2) = [ whd(1), -whd(2), -whd(3)]';
verts(1:3,3) = [-whd(1),  whd(2), -whd(3)]';
verts(1:3,4) = [ whd(1),  whd(2), -whd(3)]';
verts(1:3,5) = [-whd(1), -whd(2),  whd(3)]';
verts(1:3,6) = [ whd(1), -whd(2),  whd(3)]';
verts(1:3,7) = [-whd(1),  whd(2),  whd(3)]';
verts(1:3,8) = [ whd(1),  whd(2),  whd(3)]';
verts = E*verts;
V = verts(1:3,:)';
F = [
	1 2 4 3
	8 7 5 6
	1 2 6 5
	8 7 3 4
	1 3 7 5
	8 6 2 4
	];
end
