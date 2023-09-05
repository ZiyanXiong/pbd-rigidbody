%% Test with maximal coordinates

density = 1;
l = 5;
w = 1;
sides = [l w w];

E = eye(4);
%E(1:3,1:3) = se3.aaToMat([1 1 1],pi/4);
%E(1:3,4) = [0.5*l 0 2*l]';

%phi = [3 -4 5 -3 2 -5]';
phi = [2 2 2 0 0 0]';

E0 = E;
phi0 = phi;

clf
hold on;
axis equal;
axis(5*[-1 1 -1 1 -1 1]);
se3.drawCuboid(E,sides);
view(3);
grid on;
ax = gca;
ax.Clipping = 'off';

m = se3.inertiaCuboid(sides,density);
M = diag(m);

tEnd = 10;
h = 1e-3;
drawHz = 30;
nsteps = ceil(tEnd/h);

%%
E = E0;
phi = phi0;

t = 0;
for k = 1 : nsteps
	f = se3.ad(phi)'*M*phi;
	phi = M\(M*phi + h*f);
	E = E*se3.exp(h*phi);
	t = t + h;
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		cla;
		hold on;
		se3.drawCuboid(E,sides);
		R = E(1:3,1:3);
		p = E(1:3,4);
		w = R*phi(1:3);
		I = R*M(1:3,1:3)*R';
		l = I*w;
		ppw = [p,p+w];
		ppl = [p,p+l];
		plot3(ppw(1,:),ppw(2,:),ppw(3,:),'r');
		plot3(ppl(1,:),ppl(2,:),ppl(3,:),'g');
		title(sprintf('%f',t));
		drawnow;
	end
end

%% Test with quaternions
E = E0;
phi = phi0;

R = E(1:3,1:3);
p = E(1:3,4);
q = se3.matToQ(R);
w = R*phi(1:3); % angular velocity in world coords (omega)
v = R*phi(4:6); % linear velocity in world coords

t = 0;
for k = 1 : nsteps
	% Unconstrained step
	p0 = p;
	p = p + h*v;
	q0 = q;
	I = R*M(1:3,1:3)*R';
	wqvv = h*(I\(-cross(w,I*w)));
	w = w + wqvv;
	q = q + h*0.5*se3.qMul([w;0],q);
	q = q/norm(q);
	% Compute velocities
	v = (p - p0)/h;
	dq = se3.qMulInv(q,q0);
	w = 2*dq(1:3)/h;
	if dq(4) < 0, w = -w; end
	% Update transform
	R = se3.qToMat(q);
	t = t + h;
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		cla;
		hold on;
		E(1:3,1:3) = R;
		E(1:3,4) = p;
		se3.drawCuboid(E,sides);
		l = I*w;
		ppw = [p,p+w];
		ppl = [p,p+l];
		plot3(ppw(1,:),ppw(2,:),ppw(3,:),'r');
		plot3(ppl(1,:),ppl(2,:),ppl(3,:),'g');
		title(sprintf('%f',t));
		drawnow;
	end
end

%% Another test with quaternions
E = E0;
phi = phi0;

R = E(1:3,1:3);
p = E(1:3,4);
q = se3.matToQ(R);
pdot = R*phi(4:6);
w = phi(1:3);
W = [
	    0  w(3) -w(2)  w(1)
	-w(3)     0  w(1)  w(2)
	 w(2) -w(1)     0  w(3)
	-w(1) -w(2) -w(3)     0
	];
qdot = 0.5*W*q;
I = m(1:3);

t = 0;
for k = 1 : nsteps
	% Unconstrained step
	p0 = p;
	p = p + h*pdot;
	q0 = q;
	Q = [
		 q(4)  q(3) -q(2) -q(1)
		-q(3)  q(4)  q(1) -q(2)
		 q(2) -q(1)  q(4) -q(3)
		];
	w = 2*Q*qdot; % angular velocity in body coords
	tau = se3.cross(I.*w,w); % torque in body coords
	w = w + h*(I.\tau);
	W = [
		    0  w(3) -w(2)  w(1)
		-w(3)     0  w(1)  w(2)
		 w(2) -w(1)     0  w(3)
		-w(1) -w(2) -w(3)     0
	];
	qdot = 0.5*W*q;
	q = q + h*qdot;
	q = q/norm(q);
	% Compute velocities
	pdot = (p - p0)/h;
	qdot = (q - q0)/h;
	t = t + h;
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		cla;
		hold on;
		R = se3.qToMat(q);
		E(1:3,1:3) = R;
		E(1:3,4) = p;
		se3.drawCuboid(E,sides);
		l = I.*w;
		ppw = [p,p+R*w];
		ppl = [p,p+R*l];
		plot3(ppw(1,:),ppw(2,:),ppw(3,:),'r');
		plot3(ppl(1,:),ppl(2,:),ppl(3,:),'g');
		title(sprintf('%f',t));
		drawnow;
	end
end

%% Test with affine
E = E0;
phi = phi0;

R = E(1:3,1:3);
p = E(1:3,4);
r = reshape(R',9,1);

Edot = E*se3.brac(phi);
rdot = reshape(Edot(1:3,1:3)',9,1);
pdot = Edot(1:3,4);

R2A = 0.5*[
	-1  1  1
	 1 -1  1
	 1  1 -1
	];
I = R2A*diag(M(1:3,1:3));
W = repmat(1./I,3,1);

iters = 1;
t = 0;
for k = 1 : nsteps
	% Unconstrained step
	p0 = p;
	p = p + h*pdot;
	r0 = r;
	r = r + h*rdot;
	% Apply affine constraint
	for i = 1 : iters
		for j = 1 : 6
			[C,dC] = conAffine(r,j);
			WdC = W.*dC';
			numerator = -C;
			denominator = dC*WdC;
			dlambda = numerator/denominator;
			r = r + dlambda*WdC;
		end
	end
	% Compute velocities
	pdot = (p - p0)/h;
	rdot = (r - r0)/h;
	t = t + h;
	if drawHz > 0 && (floor(t*drawHz) > floor((t-h)*drawHz))
		cla;
		hold on;
		R = reshape(r,3,3)';
		E(1:3,1:3) = R;
		E(1:3,4) = p;
		se3.drawCuboid(E,sides);
		Rdot = reshape(rdot,3,3)';
		%phi = se3.unbrac(R'*Rdot);
		phibrac = R'*Rdot;
		phi = se3.unbrac(0.5*(phibrac - phibrac'));
		w = R*phi(1:3);
		I = R*M(1:3,1:3)*R';
		l = I*w;
		ppw = [p,p+w];
		ppl = [p,p+l];
		plot3(ppw(1,:),ppw(2,:),ppw(3,:),'r');
		plot3(ppl(1,:),ppl(2,:),ppl(3,:),'g');
		title(sprintf('%f',t));
		drawnow;
	end
end

%%
function [C,dC] = conAffine(r,i)
A11 = r(1);
A12 = r(2);
A13 = r(3);
A21 = r(4);
A22 = r(5);
A23 = r(6);
A31 = r(7);
A32 = r(8);
A33 = r(9);
switch i
	case 1
		C = A11^2 + A21^2 + A31^2 - 1;
		dC = [2*A11, 0, 0, 2*A21, 0, 0, 2*A31, 0, 0];
	case 2
		C = A12^2 + A22^2 + A32^2 - 1;
		dC = [0, 2*A12, 0, 0, 2*A22, 0, 0, 2*A32, 0];
	case 3
		C = A13^2 + A23^2 + A33^2 - 1;
		dC = [0, 0, 2*A13, 0, 0, 2*A23, 0, 0, 2*A33];
	case 4
		C = A11*A12 + A21*A22 + A31*A32;
		dC = [A12, A11, 0, A22, A21, 0, A32, A31, 0];
	case 5
		C = A12*A13 + A22*A23 + A32*A33;
		dC = [0, A13, A12, 0, A23, A22, 0, A33, A32];
	case 6
		C = A11*A13 + A21*A23 + A31*A33;
		dC = [A13, 0, A11, A23, 0, A21, A33, 0, A31];
end
end
