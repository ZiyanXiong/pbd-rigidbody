scene = Scene();
density = 1.0;
sides = [2 1 0.5];
scene.ground.E = eye(4);
scene.tEnd = 3;
scene.bodies{1} = BodyCuboid(density,sides); %#ok<*SAGROW>
scene.bodies{1}.v = [0 0 0]';
scene.bodies{1}.w = [0 0 0]';
scene.bodies{1}.x = [0 0 0.25]';
scene.bodies{1}.q = quaternion([1 0.0 0.0 0.0]);
scene.bodies{1}.updateE();
scene.bodies{1}.collide = true;
scene.bodies{1}.name = 'body1';
scene.bodies{2} = BodyCuboid(density,sides); %#ok<*SAGROW>
scene.bodies{2}.v = [0 0 0]';
scene.bodies{2}.w = [0 0 0]';
scene.bodies{2}.x = [0.25 0 0.75]';
scene.bodies{2}.q = quaternion([1 0.0 0.0 0.0]);
scene.bodies{2}.updateE();
scene.bodies{2}.collide = true;
scene.bodies{2}.name = 'body2';
scene.constraints{1} = ConstraintGroundContact({scene.bodies{1}}, scene.ground.E);
%scene.constraints{2} = ConstraintSphericalJoint({scene.bodies{1}}, [1.0 0.0 0]');
%scene.constraints{2} = ConstraintSphericalJoint({scene.bodies{2}}, [1.0 0.0 0]');
scene.constraints{2} = ConstraintGroundContact({scene.bodies{2}}, scene.ground.E);
scene.constraints{3} = ConstraintBodiesContact({scene.bodies{1}, scene.bodies{2}});

scene.init();
grav = [0 0 -9.8]';
%grav = [0 0 0]';

nsteps = scene.nsteps;
for i = 0 : nsteps - 1
    h = scene.h / scene.sub_steps;
    for j = 0 : scene.sub_steps -1
        %scene.collide();
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateWithoutConstraints(grav, h);
        end
        scene.solvePositions();
        %scene.bodies{1}.v
        for k = 1 : length(scene.bodies)
            scene.bodies{k}.updateAfterSolve(h);
        end
        %scene.bodies{1}.v
        %draw(scene);
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
for i = 1:length(scene.bodies)
    E = scene.bodies{i}.E_wi;
    sides = scene.bodies{i}.sides;
    se3.drawAxis(E);
    se3.drawCuboid(E,sides);
end
    title(sprintf('%f',scene.t));
    drawnow
end

