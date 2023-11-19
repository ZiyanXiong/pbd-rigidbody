% /******************************************************************************
% 
%   This source code accompanies the Journal of Graphics Tools paper:
% 
%   "Fast Ray-Axis Aligned Bounding Box Overlap Tests With Pluecker Coordinates" by
%   Jeffrey Mahovsky and Brian Wyvill
%   Department of Computer Science, University of Calgary
% 
%   This source code is public domain, but please mention us if you use it.
% 
%  ******************************************************************************/
function [hit,t] = smits_mul(rx, rv, bmin, bmax)
% rx: ray origin (3x1)
% rv: ray direction (3x1)
% bmin: box min (3x1)
% bmax: box max (3x1)

r = make_ray(rx(1), rx(2), rx(3), rv(1), rv(2), rv(3));
b = make_aabox(bmin(1), bmin(2), bmin(3), bmax(1), bmax(2), bmax(3));

tnear = -1e6;
tfar = 1e6;

hit = false;
t = 0;

% multiply by the inverse instead of dividing
t1 = (b.x0 - r.x) * r.ii;
t2 = (b.x1 - r.x) * r.ii;

if t1 > t2
	temp = t1;
	t1 = t2;
	t2 = temp;
end
if t1 > tnear
	tnear = t1;
end
if t2 < tfar
	tfar = t2;
end
if tnear > tfar
	return;
end
if tfar < 0.0
	return;
end

t1 = (b.y0 - r.y) * r.ij;
t2 = (b.y1 - r.y) * r.ij;

if t1 > t2
	temp = t1;
	t1 = t2;
	t2 = temp;
end
if t1 > tnear
	tnear = t1;
end
if t2 < tfar
	tfar = t2;
end
if tnear > tfar
	return;
end
if tfar < 0.0
	return
end

t1 = (b.z0 - r.z) * r.ik;
t2 = (b.z1 - r.z) * r.ik;

if t1 > t2
	temp = t1;
	t1 = t2;
	t2 = temp;
end
if t1 > tnear
	tnear = t1;
end
if t2 < tfar
	tfar = t2;
end
if tnear > tfar
	return
end
if tfar < 0.0
	return
end

t = tnear;
hit = true;

end

%%
function r = make_ray(x, y, z, i, j, k)
r.x = x;
r.y = y;
r.z = z;
r.i = i;
r.j = j;
r.k = k;
r.ii = 1.0 / i;
r.ij = 1.0 / j;
r.ik = 1.0 / k;
end

%%
function a = make_aabox(x0, y0, z0, x1, y1, z1)
if x0 > x1
	a.x0 = x1;
	a.x1 = x0;
else
	a.x0 = x0;
	a.x1 = x1;
end
if y0 > y1
	a.y0 = y1;
	a.y1 = y0;
else
	a.y0 = y0;
	a.y1 = y1;
end
if z0 > z1
	a.z0 = z1;
	a.z1 = z0;
else
	a.z0 = z0;
	a.z1 = z1;
end
end
