% https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#From_a_quaternion_to_an_orthogonal_matrix
clear;
syms wx wy wz w0 real
w00 = w0*w0;
w0x = w0*wx;
w0y = w0*wy;
w0z = w0*wz;
wxx = wx*wx;
wxy = wx*wy;
wzx = wz*wx;
wyy = wy*wy;
wyz = wy*wz;
wzz = wz*wz;
R = [
	w00 + wxx - wyy - wzz, 2*(wxy - w0z), 2*(wzx + w0y)
	2*(wxy + w0z), w00 - wxx + wyy - wzz, 2*(wyz - w0x)
	2*(wzx - w0y), 2*(wyz + w0x), w00 - wxx - wyy + wzz
	];
r = reshape(R',9,1);
Jaw = jacobian(r,[wx wy wz w0]);
syms mx my mz real
I = Jaw'*diag([mx my mz mx my mz mx my mz])*Jaw;
syms dwx dwy dwz dw0 real
dJaw = 2*[
	 dwx, -dwy, -dwz,  dw0
	 dwy,  dwx, -dw0, -dwz
	 dwz,  dw0,  dwx,  dwy
	 dwy,  dwx,  dw0,  dwz
	-dwx,  dwy, -dwz,  dw0
	-dw0,  dwz,  dwy, -dwx
	 dwz, -dw0,  dwx, -dwy
	 dw0,  dwz,  dwy,  dwx
	-dwx, -dwy,  dwz,  dw0
	];
fqvv = -Jaw'*diag([mx my mz mx my mz mx my mz])*dJaw*[dwx dwy dwz dw0]';

maxyz = mx + my + mz;
maxy = mx + my;
mayz = my + mz;
mazx = mz + mx;
msxz = mx - mz;
msyx = my - mx;
mszy = mz - my;
M(1,:) = [mayz*w00 + maxyz*wxx + maxy*wyy + mazx*wzz, msyx*w0z + mz*wxy, msxz*w0y + my*wzx, mx*w0x + mszy*wyz];
M(2,:) = [M(1,2), mazx*w00 + maxy*wxx + maxyz*wyy + mayz*wzz, mszy*w0x + mx*wyz, my*w0y + msxz*wzx];
M(3,:) = [M(1,3), M(2,3), maxy*w00 + mazx*wxx + maxyz*wzz + mayz*wyy, mz*w0z + msyx*wxy];
M(4,:) = [M(1,4), M(2,4), M(3,4), maxyz*w00 + mazx*wyy + mayz*wxx + maxy*wzz];
M = 4*M;
simplify(M - I)

dw0x = dw0*dwx;
dw0y = dw0*dwy;
dw0z = dw0*dwz;
dw00 = dw0*dw0;
dwxx = dwx*dwx;
dwxy = dwx*dwy;
dwyy = dwy*dwy;
dwyz = dwy*dwz;
dwzx = dwx*dwz;
dwzz = dwz*dwz;
fxx = -dw00 - dwxx + dwyy + dwzz;
fxy =  dw00 - dwxx + dwyy - dwzz;
fxz =  dw00 - dwxx - dwyy + dwzz;
fyx =  dw00 + dwxx - dwyy - dwzz;
fyy = -dw00 + dwxx - dwyy + dwzz;
fyz =  dw00 - dwxx - dwyy + dwzz;
fzx =  dw00 + dwxx - dwyy - dwzz;
fzy =  dw00 - dwxx + dwyy - dwzz;
fzz = -dw00 + dwxx + dwyy - dwzz;
f0x = -dw00 - dwxx + dwyy + dwzz;
f0y = -dw00 + dwxx - dwyy + dwzz;
f0z = -dw00 + dwxx + dwyy - dwzz;
mxdwxya0z = mx*(dwxy + dw0z);
mydwxys0z = my*(dwxy - dw0z);
mzdwzxa0y = mz*(dwzx + dw0y);
mxdwzxs0y = mx*(dwzx - dw0y);
mydwyza0x = my*(dwyz + dw0x);
mzdwyzs0x = mz*(dwyz - dw0x);
fxax = mydwyza0x + mzdwyzs0x;
fxsx = mydwyza0x - mzdwyzs0x;
fyay = mzdwzxa0y + mxdwzxs0y;
fysy = mzdwzxa0y - mxdwzxs0y;
fzaz = mxdwxya0z + mydwxys0z;
fzsz = mxdwxya0z - mydwxys0z;
f(1,1) = wx*(mx*fxx + my*fxy + mz*fxz) - 2*(          wy*fzaz + wz*fyay + w0*fxsx);
f(2,1) = wy*(mx*fyx + my*fyy + mz*fyz) - 2*(wx*fzaz           + wz*fxax + w0*fysy);
f(3,1) = wz*(mx*fzx + my*fzy + mz*fzz) - 2*(wx*fyay + wy*fxax           + w0*fzsz);
f(4,1) = w0*(mx*f0x + my*f0y + mz*f0z) - 2*(wx*fxsx + wy*fysy + wz*fzsz);
f = 4*f;
simplify(f - fqvv)

%% Transform between quaternion derivative and angular velocity
% https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
% http://www.iri.upc.edu/people/jsola/JoanSola/objectes/notes/kinematics.pdf
clear;
q = randn(4,1); % quaternion, rotation wrt *world space*
q = q/norm(q);
if q(4) < 0, q = -q; end
w = randn(3,1); % angular velocity in *body space*
qdot = 0.5*[
	           + w(3)*q(2) - w(2)*q(3) + w(1)*q(4) 
	-w(3)*q(1)             + w(1)*q(3) + w(2)*q(4) 
	 w(2)*q(1) - w(1)*q(2)             + w(3)*q(4)
	-w(1)*q(1) - w(2)*q(2) - w(3)*q(3)];
% Finite difference check
h = 1e-6;
R0 = se3.qToMat(q);
R1 = R0*se3.exp(h*w(1:3));
q0 = se3.matToQ(R0);
q1 = se3.matToQ(R1);
if q0(4) < 0, q0 = -q0; end
if q1(4) < 0, q1 = -q1; end
qdot_ = (q1 - q0)/h;
fprintf('%e\n',norm(qdot_ - qdot));

% Using a matrix: qdot = 0.5*Omega(w)*q
%            [    0  w(3) -w(2)  w(1)]
% Omega(w) = [-w(3)     0  w(1)  w(2)]
%            [ w(2) -w(1)     0  w(3)]
%            [-w(1) -w(2) -w(3)     0]
%O = [-se3.brac(w(1:3)), w(1:3); -w(1:3)', 0];
O = [
	    0  w(3) -w(2)  w(1)
	-w(3)     0  w(1)  w(2)
	 w(2) -w(1)     0  w(3)
	-w(1) -w(2) -w(3)     0
	];
fprintf('%e\n',norm(qdot - 0.5*O*q));

% Inverse relationship: get phi from wdot
P = [
	 q(4) -q(3)  q(2)
	 q(3)  q(4) -q(1)
	-q(2)  q(1)  q(4)
	-q(1) -q(2) -q(3)
	];
% qdot = 0.5*P*w
% 2*P'*qdot = P'*P*w
% w = 2*P'*qdot, since P'*P = I
fprintf('%e\n',norm(w - 2*P'*qdot));

% Non-matrix version
w1 = 2*[
	 q(4)*qdot(1) + q(3)*qdot(2) - q(2)*qdot(3) - q(1)*qdot(4)
	-q(3)*qdot(1) + q(4)*qdot(2) + q(1)*qdot(3) - q(2)*qdot(4)
	 q(2)*qdot(1) - q(1)*qdot(2) + q(4)*qdot(3) - q(3)*qdot(4)
	];
fprintf('%e\n',norm(w - w1));

% If we want to apply a torque in body space, we need the Jacobian
% transpose.
%    w = J*qdot
%    f = J'*t
% where t is the body-space torque, and f is the generalized force acting
% on the quaternions.
J = 2*P';
