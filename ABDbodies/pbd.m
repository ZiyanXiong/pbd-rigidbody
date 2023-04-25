scene = Scene();
density = 1.0;
sides = [2 1 0.5];
scene.ground.E = eye(4);
scene.tEnd = 5;
scene.bodies{1} = BodyCuboid(density,sides); %#ok<*SAGROW>
%E = eye(4);
E = se3.randE();
%E(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
E(1:3,4) = [0 0 4]';
scene.bodies{1}.setBodyTransform(E(1:3,4), E(1:3,1:3));
scene.bodies{1}.phi = [10 1 0 0 0 0]'; % rigid velocity in local space
scene.bodies{1}.updateQdot();
scene.bodies{1}.collide = true;
scene.constraints{1} = ConstraintGroundContact({scene.bodies{1}}, scene.ground.E);
scene.init();
grav = [0 0 -9.8]';
%grav = [0 0 0]';

f = zeros(12,1);
f(1:3) = scene.bodies{1}.M(1:3).*grav;

nsteps = scene.nsteps;
for i = 0 : nsteps - 1
    h = scene.h / scene.sub_steps;
    for j = 0 : scene.sub_steps -1
        %scene.collide();
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateWithoutConstraints(f, h);
        end
        scene.solvePositions();
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateAfterSolve(h);
        end
        %scene.bodies{1}.w
        draw(scene);
    end
    %draw(scene);
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

