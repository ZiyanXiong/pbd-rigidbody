rng(42);
%[w,x,y,z] = parts(randrot);
%q1 = [x, y, z, w]';
%[w,x,y,z] = parts(randrot);
%q2 = [x, y, z, w]';

q1 = [ 0, -0.3826834, 0, 0.9238795 ]';
q2 = [ 0.2209424, 0.2209424, 0.2209424, 0.9238795 ]';


axis = [0 1 0]';
b1 = [1 0 0]';
b2 = [0 0 1]';
% cf = [se3.qRot(q1, axis); se3.qRot(q1, b1); se3.qRot(q1, b2)];
dtheta1 = asin(cross(se3.qRot(q2, axis), se3.qRot(q1, axis)));
% dtheta2 = asin(cross(se3.qRot(q2, b1), se3.qRot(q1, b1)));
% dtheta3 = cross(se3.qRot(q2, b2), se3.qRot(q1, b2));

q_mid_1 = se3.qMul(se3.deltaThetaToDq(dtheta1), q2);
dtheta2 = asin(cross(se3.qRot(q_mid_1, b1), se3.qRot(q1, b1)));
q_mid_2 = se3.qMul(se3.deltaThetaToDq(dtheta2), q_mid_1);
q_end = se3.qMul(se3.qMul(se3.deltaThetaToDq(dtheta2), se3.deltaThetaToDq(dtheta1)), q2);
theta_end = se3.dqToDeltaTheta(se3.qMul(se3.deltaThetaToDq(dtheta2), se3.deltaThetaToDq(dtheta1)));
dq = se3.qMulInv(q1, q_end);

vo1 = getOrthogonal(axis);
vo2 = getOrthogonal(se3.qRot(q2, axis));
vo3 = se3.qRotInv(q2,vo2);

function vo = getOrthogonal(v)
    g = sign(v(3));
    if(g==0)
        g = 1;
    end
    h = v(3) + g;
    vo = [g - v(1)*v(1)/h, -v(1)*v(2)/h, -v(1)]';
end

