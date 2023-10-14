function model = createTestModels(modelID)

model = apbd.Model();

switch(modelID)
	case 0
		model.name = 'Rigid Body';
		model.plotH = true;
		model.tEnd = 1;
		model.h = 1/30;
		model.substeps = 10;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = 0*[0 0 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density); %#ok<*SAGROW>
		E = eye(4);
		R = se3.aaToMat([1 1 1],pi/4);
		E(1:3,1:3) = R;
		E(1:3,4) = [0 0 5]';
		model.bodies{end}.setInitTransform(E);
		model.bodies{end}.setInitVelocity([R'*[3 -4 5]'; R'*[0 0 5]']);
		%model.bodies{end}.setInitVelocity([0 0 0 0 0 1]');
	case 1
		model.name = 'Spherical Joint';
		model.plotH = true;
		model.tEnd = 1;
		model.h = 1/30;
		model.substeps = 20;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = [-10 10 -10 10 0 10];
		%model.view = [0 0];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		E = eye(4);
		E(1:3,4) = [0.5*l 0 2*l]';
		model.bodies{end}.setInitTransform(E);
		model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
		model.constraints{end+1} = apbd.ConJointSph({model.bodies{end}},{[-l/2 0 0]'}); %#ok<*CCAT1>
		model.constraints{end}.radius = w/2;

		E = eye(4);
		E(1:3,4) = [1.5*l 0 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		model.constraints{end+1} = apbd.ConJointSph({model.bodies{end-1},model.bodies{end}},{[l/2 0 0]',[-l/2 0 0]'});
		model.constraints{end}.radius = w/2;
	case 2
		model.name = 'Fixed Joint';
		model.plotH = true;
		model.tEnd = 1;
		model.h = 1/30;
		model.substeps = 20;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = [-10 10 -10 10 0 10];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		E = eye(4);
		E(1:3,4) = [0.5*l 0 2*l]';
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end}},{E1});
		E1(1:3,4) = [l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end}},{E1});

		E = eye(4);
		E(1:3,4) = [1.5*l 0 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,4) = [l/2 0 0]';
		E2 = eye(4);
		E2(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end-1},model.bodies{end}},{E1,E2});
	case 3
		model.name = 'Revolute Joint';
		model.plotH = true;
		model.tEnd = 1.0;
		model.h = 1/30;
		model.substeps = 30;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = [0 -980 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = [-10 10 -10 10 0 10];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		E = eye(4);
		E(1:3,4) = [0.5*l 0 2*l]';
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E1(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointRev({model.bodies{end}},{E1});
		model.constraints{end}.radius = 0.25;
		model.constraints{end}.height = 1.5;

		E = eye(4);
		%E(1:3,1:3) = se3.aaToMat([1 0 0],pi/6);
		E(1:3,4) = [1.5*l 0 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E1(1:3,4) = [l/2 0 0]';
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E2(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointRev({model.bodies{end-1},model.bodies{end}},{E1,E2});
		model.constraints{end}.radius = 0.25;
		model.constraints{end}.height = 1.5;
	case 4
		model.name = 'Prismatic Joint';
		model.plotH = true;
		model.tEnd = 0.55;
		model.h = 1e-3;
		model.substeps = 10;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = [10 10 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = [-10 10 -10 10 0 10];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		E = eye(4);
		E(1:3,4) = [0.5*l 0 2*l]';
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		E1(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end}},{E1});
		%model.constraints{end}.width = 0.5;
		%model.constraints{end}.height = 2.0;
		E1(1:3,4) = [l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end}},{E1});

		E = eye(4);
		E(1:3,4) = [1.5*l 0 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		E1(1:3,4) = [l/2 0 0]';
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		E2(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointPri({model.bodies{end-1},model.bodies{end}},{E1,E2});
		model.constraints{end}.width = 0.5;
		model.constraints{end}.height = 2.0;

		E = eye(4);
		E(1:3,1:3) = se3.aaToMat([0 0 1],pi/2);
		E(1:3,4) = [1.5*l 0.5*l 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([0 0 1],-pi/2);
		E2(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end-1},model.bodies{end}},{E1,E2});
	case 5
		model.name = 'Cylindrical Joint';
		model.plotH = true;
		model.tEnd = 0.55;
		model.h = 1/30;
		model.substeps = 30;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = [10 10 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = [-10 10 -10 10 0 10];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		E = eye(4);
		E(1:3,4) = [0.5*l 0 2*l]';
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		E1(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end}},{E1});
		%model.constraints{end}.radius = 0.25;
		%model.constraints{end}.height = 2.0;
		E1(1:3,4) = [l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end}},{E1});

		E = eye(4);
		E(1:3,4) = [1.5*l 0 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		E1(1:3,4) = [l/2 0 0]';
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		E2(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointCyl({model.bodies{end-1},model.bodies{end}},{E1,E2});
		model.constraints{end}.radius = 0.25;
		model.constraints{end}.height = 2.0;

		E = eye(4);
		E(1:3,1:3) = se3.aaToMat([0 0 1],pi/2);
		E(1:3,4) = [1.5*l 0.5*l 2*l]';
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
		E1 = eye(4);
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([0 0 1],-pi/2);
		E2(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointFix({model.bodies{end-1},model.bodies{end}},{E1,E2});
	case 6
		model.name = 'Ground Collision';
		model.plotH = false;
		model.tEnd = 1.0;
		model.h = 1e-3;
		model.substeps = 1;
		model.iters = 2;
		density = 1.0;
		l = 1;
		w = 1;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		angle = 20*pi/180;
		R = se3.aaToMat([0 1 0],angle);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 10;
		model.axis = 5*[-1 1 -1 1 0 1];
		model.drawHz = 1000;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = 1.0*(sin(angle)/cos(angle));
		E = eye(4);
		E(1:3,1:3) = R;
		E(1:3,4) = R*[0 0 w/2]';
		%E(1:3,4) = [0 0 5]';
		model.bodies{end}.setInitTransform(E);
		%model.bodies{end}.setInitVelocity([0 0 0 0 0 10]');
	case 7
		model.name = 'Rigid Collisions';
		model.plotH = false;
		model.tEnd = 1;
		model.h = 1e-2;
		model.substeps = 1;
		model.iters = 200;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.9;

		model.ground.size = 10;
		model.axis = 5*[-1 1 -1 1 0 1];
		model.drawHz = 10000;

		model.view = [0 0];

		n = 2;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
			E = eye(4);
			x = 0.1*i;
			y = 0;
			z = (i-0.5)*w*0.99;
			E(1:3,4) = [x y z]';
			model.bodies{end}.setInitTransform(E);
            model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
		end
	case 8
		model.name = '2D Rigid Body';
		model.plotH = false;
		model.tEnd = 10.5;
		model.h = 1e-2;
		model.substeps = 1;
		model.iters = 10;
		density = 1.0;
		l = 1;
		w = 1;
		sides = [l w w];
		model.grav = [0 -980 0]';
		model.ground.E = eye(4);

		R = se3.aaToMat([1 0 0],-pi/2);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 10;
		model.axis = 5*[-1 1 -0 1 -1 1];
		model.view = 2;
		model.drawHz = 10000;
		n = 1;
		for i = 1 : n
		    model.bodies{end+1} = apbd.BodyRigid2d(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = 0.9;
            if(i ~= 1)
                %model.bodies{end}.color = model.bodies{end - 1}.color;
            end
		    E = eye(4);
			x = 0.05*i;
			y = (i-0.5)*w*0.99;
			z = 0;
			E(1:3,4) = [x y z]';
            %E(1:3,4) = [0 -0.5 0]';
		    model.bodies{end}.setInitTransform(E);
		    model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
        end
	case 9
		model.name = 'Affine Body';
		model.plotH = true;
		model.tEnd = 1;
		model.h = 1/30;
		model.substeps = 10;
		model.iters = 1;
		density = 1.0;
		l = 5;
		w = 1;
		sides = [l w w];
		model.grav = 0*[0 0 -980]';
		model.ground.E = eye(4);

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 30;

		model.bodies{end+1} = apbd.BodyAffine(apbd.ShapeCuboid(sides),density);
		model.constraints{end+1} = apbd.ConAffine(model.bodies{end});
		E = eye(4);
		R = se3.aaToMat([1 1 1],pi/4);
		E(1:3,1:3) = R;
		E(1:3,4) = [0 0 5]';
		model.bodies{end}.setInitTransform(E);
		model.bodies{end}.setInitVelocity([R'*[3 -4 5]'; R'*[0 0 5]']);
		%model.bodies{end}.setInitVelocity([2 2 2 0 0 0]');
	case 10
		model.name = 'Affine Ground Collision';
		model.plotH = false;
		model.tEnd = 1.0;
		model.h = 1e-3;
		model.substeps = 1;
		model.iters = 5; % iters and substeps are not equivalent!
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		angle = 10*pi/180;
		R = se3.aaToMat([0 1 0],angle);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 10;
		model.axis = 5*[-1 1 -1 1 0 1];
		model.drawHz = 30;

		if true
			% Affine
			model.bodies{end+1} = apbd.BodyAffine(apbd.ShapeCuboid(sides),density);
			model.constraints{end+1} = apbd.ConAffine(model.bodies{end});
		else
			% Rigid for testing
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density); %#ok<UNRCH>
		end
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = 0.8;%(sin(angle)/cos(angle));
		E = eye(4);
		%E(1:3,1:3) = R;
		%E(1:3,4) = R*[0 0 w/2]';
		E(1:3,4) = [0 0 5]';
		model.bodies{end}.setInitTransform(E);
		%model.bodies{end}.setInitVelocity([0 0 0 0 0 10]');

        case 11
		model.name = 'Affine Body Ground';
		model.plotH = false;
		model.tEnd = 10.5;
		model.h = 1e-2;
		model.substeps = 1;
		model.iters = 400;
		density = 1.0;
		l = 1;
		w = 1;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		model.ground.size = 10;

        %model.axis = 5*[-1 1 -0 1 -1 1];
		model.axis = 5*[-1 1 -1 1 0 1];
		model.drawHz = 10000;
		model.view = [0 0];
		n = 1;
		for i = 1 : n
    		model.bodies{end+1} = apbd.BodyAffine(apbd.ShapeCuboid(sides),density);
		    %model.constraints{end+1} = apbd.ConAffine(model.bodies{end});
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = 0.9;
		    E = eye(4);
			x = 0.1*i;
			y = 0;
			z = (i-0.5)*w*0.99;
			E(1:3,4) = [x y z]';
		    model.bodies{end}.setInitTransform(E);
		    model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
        end
end

end
