scene = Scene();
density = 1.0;
sides = [2 1 0.5];
scene.tEnd = 5;
scene.bodies{1} = BodyCuboid(density,sides); %#ok<*SAGROW>
scene.bodies{1}.v = [0 0 0]';
scene.bodies{1}.w = [0 0 0]';
scene.bodies{1}.x = [0 0 2]';
scene.bodies{1}.q = quaternion([1 0 0.2 0]);
scene.bodies{1}.E_wi(1:3,1:3) = quat2rotm(scene.bodies{1}.q);
scene.bodies{1}.E_wi(1:3,4) = scene.bodies{1}.x;
scene.bodies{1}.collide = true;
scene.ground.E = eye(4);
scene.init();
grav = [0 0 -9.8]';
%grav = [0 0 0]';

nsteps = scene.nsteps;
for i = 0 : nsteps - 1
    h = scene.h / scene.sub_steps;
    for j = 0 : scene.sub_steps -1
        scene.collide();
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateWithoutConstraints(grav, h);
        end
        scene.solvePositions();
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateAfterSolve(h);
        end
    end
    draw(scene);
    scene.t = scene.t + scene.h;
end

%%
function draw(scene)
if scene.t == 0
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
E = scene.bodies{1}.E_wi;
sides = scene.bodies{1}.sides;
se3.drawAxis(E);
se3.drawCuboid(E,sides);
title(sprintf('%f',scene.t));
drawnow
end

