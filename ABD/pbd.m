clear;
rng(2);

scene = Scene();
density = 1.0;
sides = [2 1 0.5];
scene.ground.E = eye(4);
scene.tEnd = 5;
scene.bodies{1} = BodyCuboid(density,sides); %#ok<*SAGROW>
scene.bodies{2} = BodyCuboid(density,sides); %#ok<*SAGROW>
%scene.bodies{3} = BodyCuboid(density,sides); %#ok<*SAGROW>
%scene.bodies{4} = BodyCuboid(density,sides); %#ok<*SAGROW>
E = eye(4);
%E = se3.randE();
%E(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
E(1:3,4) = [0 0 0.25]';
scene.bodies{1}.setBodyTransform(E(1:3,4), E(1:3,1:3));
scene.bodies{1}.phi = [0 0 0 0 0 0]'; % rigid velocity in local space
scene.bodies{1}.updateQdot();
scene.bodies{1}.collide = true;

E = eye(4);
E(1:3,4) = [0.25 0 0.75]';
scene.bodies{2}.setBodyTransform(E(1:3,4), E(1:3,1:3));
scene.bodies{2}.phi = [0 0 0 0 0 0]'; % rigid velocity in local space
scene.bodies{2}.updateQdot();
scene.bodies{2}.collide = true;

%{
E = eye(4);
E(1:3,4) = [0.1 0 1.25]';
scene.bodies{3}.setBodyTransform(E(1:3,4), E(1:3,1:3));
scene.bodies{3}.phi = [0 0 0 0 0 0]'; % rigid velocity in local space
scene.bodies{3}.updateQdot();
scene.bodies{3}.collide = true;

E = eye(4);
E(1:3,4) = [0.15 0 1.75]';
scene.bodies{4}.setBodyTransform(E(1:3,4), E(1:3,1:3));
scene.bodies{4}.phi = [0 0 0 0 0 0]'; % rigid velocity in local space
scene.bodies{4}.updateQdot();
scene.bodies{4}.collide = true;
%}


%scene.constraints{7} = ConstraintVolume({scene.bodies{4}});
%scene.constraints{8} = ConstraintOrtho({scene.bodies{4}});



scene.constraints{5} = ConstraintBodiesContact({scene.bodies{1}, scene.bodies{2}});
scene.constraints{3} = ConstraintVolume({scene.bodies{1}});
scene.constraints{4} = ConstraintOrtho({scene.bodies{1}});
scene.constraints{1} = ConstraintVolume({scene.bodies{2}});
scene.constraints{2} = ConstraintOrtho({scene.bodies{2}});

scene.constraints{6} = ConstraintGroundContact({scene.bodies{1}}, scene.ground.E);
scene.constraints{3} = ConstraintVolume({scene.bodies{1}});
scene.constraints{4} = ConstraintOrtho({scene.bodies{1}});

%{
scene.constraints{9} = ConstraintBodiesContact({scene.bodies{3}, scene.bodies{2}});

scene.constraints{10} = ConstraintVolume({scene.bodies{2}});
scene.constraints{11} = ConstraintOrtho({scene.bodies{2}});
scene.constraints{12} = ConstraintVolume({scene.bodies{3}});
scene.constraints{13} = ConstraintOrtho({scene.bodies{3}});
%}

%scene.constraints{11} = ConstraintBodiesContact({scene.bodies{2}, scene.bodies{3}});
%scene.constraints{12} = ConstraintBodiesContact({scene.bodies{3}, scene.bodies{4}});
%scene.constraints{9} = ConstraintGroundContact({scene.bodies{2}}, scene.ground.E);
%scene.constraints{7} = ConstraintSphericalJoint({scene.bodies{1}}, [1.0 0.5 0.0]');
%scene.constraints{7} = ConstraintSphericalJoint({scene.bodies{1}}, [-1.0 0.5 0.0]');
%scene.constraints{8} = ConstraintSphericalJoint({scene.bodies{1}}, [-1.0 -0.5 0.0]');
%scene.constraints{9} = ConstraintSphericalJoint({scene.bodies{1}}, [0.0 0.0 0.0]');
%scene.constraints{8} = ConstraintSphericalJoint({scene.bodies{1}}, [-1.0 -0.5 0.0]');
%scene.constraints{9} = ConstraintSphericalJoint({scene.bodies{1}, scene.bodies{2}}, [1.0 0.0 0.0]');

scene.sub_steps = 1;
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
for i = 1 : length(scene.bodies)
    E = scene.bodies{i}.E_wi;
    sides = scene.bodies{i}.sides;
    se3.drawAxis(E);
    se3.drawCuboid(E,sides);
end
title(sprintf('%f',scene.t));
drawnow
end

