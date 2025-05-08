% Parameters
dt = 0.01;             % Time step
steps = 100;           % Number of steps
kp = 1000;             % Proportional gain
kd = 100;              % Derivative gain

% Initial state
pos = 0.0;             % Initial position
vel = 0.0;             % Initial velocity

% Target state
target_pos = pi/4;     % Target position
target_vel = 0.0;      % Target velocity

% Storage for plotting
positions_error = zeros(1, steps);

% Simulation loop
for i = 1:steps
    %[pos,vel] = explicitDrive(pos,vel,target_pos, target_vel, kp, kd, dt);
    %[pos,vel] = implicitDrive(pos,vel,target_pos, target_vel, kp, kd, dt);
    [pos,vel] = semiImplicitDrive(pos,vel,target_pos, target_vel, kp, kd, dt);
    positions_error(i) = target_pos - pos;
end

% Time vector
time = (0:steps-1) * dt;

% Plot
figure;
plot(time, positions_error, 'b', 'LineWidth', 2);
hold on;
xlabel('Time (s)');
ylabel('Position Errors');
title('PD-Controlled Particle (PBD Style)');
grid on;

function [pos,vel] = explicitDrive(pos,vel,xt,vt,kp,kd, dt)
    lambda = dt*(kp*(xt-pos-vel*dt) + kd*(vt-vel));
    vel = vel + lambda;
    pos = pos + vel * dt;
end

function [pos,vel] = implicitDrive(pos,vel,xt,vt,kp,kd, dt)
    a = dt*(dt*kp+kd);
    b = dt*kd*vt;
    x = 1/(a+1);
    itermax = 1;
    lambda = 0;
    for i = 1:itermax
        lambda = x*b + x*dt*kp*(xt-pos) - x*a*vel + (1-x)*lambda;
        vel = vel + lambda;
    end
    pos = pos + vel * dt;
end

function [pos,vel] = semiImplicitDrive(pos,vel,xt,vt,kp,kd, dt)
    a = dt*(dt*kp+kd);
    x = 1/(a+1);
    lambda = x*dt*(kp*(xt-pos-vel*dt) + kd*(vt-vel));
    vel = vel + lambda;
    pos = pos + vel * dt;
end