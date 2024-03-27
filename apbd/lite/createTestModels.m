function model = createTestModels(modelID)

model = apbd.Model();

switch(modelID)
	case 0
		model.name = 'Rigid Collisions';
		model.plotH = false;
		model.tEnd = 1;
		model.h = 1e-2;
		model.substeps = 1;
		model.iters = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

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
			x = 0.25*i;
			y = 0;
			z = (i-0.5)*w*0.99+1;
			E(1:3,4) = [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end


    	case 1
		model.name = 'Rigid Collisions with Joints and stacks';
		model.plotH = false;
		model.tEnd = 1;
		model.h = 5e-3;
		model.substeps = 1;
		model.iters = 30;
		density = 1.0;
		w = 1;
        l = 2;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.3;

		model.ground.size = 10;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 10000;

		model.view = [0 0];
        radius = 0.6;
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
	    model.bodies{end}.collide = false;
	    model.bodies{end}.mu = mu;
		E = eye(4);
		E(1:3,4) = [-1.5*l-3*radius 0 3*l+radius]';
		model.bodies{end}.setInitTransform(E);

		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E1(1:3,4) = [-l/2 0 0]';
		model.constraints{end+1} = apbd.ConJointRev({model.bodies{end}},{E1});
		model.constraints{end}.radius = radius*2;
		model.constraints{end}.height = 1.5;

		E = eye(4);
		%E(1:3,1:3) = se3.aaToMat([1 0 0],pi/6);
		E(1:3,4) = [-0.5*l-radius 0 3*l+radius]';
		%E(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
	    model.bodies{end}.collide = false;
	    model.bodies{end}.mu = mu;
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E1(1:3,4) = [l/2+radius 0 0]';
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E2(1:3,4) = [-l/2-radius 0 0]';
		model.constraints{end+1} = apbd.ConJointRev({model.bodies{end},model.bodies{end-1}},{E2,E1});
		model.constraints{end}.radius = radius*2;
		model.constraints{end}.height = 1.5;

		E = eye(4);
		%E(1:3,1:3) = se3.aaToMat([1 0 0],pi/6);
		E(1:3,4) = [0 0 2.5*l]';
		E(1:3,1:3) = se3.aaToMat([0 1 0],pi/2);
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.setInitTransform(E);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
		E1 = eye(4);
		E1(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E1(1:3,4) = [l/2+radius 0 0]';
		E2 = eye(4);
		E2(1:3,1:3) = se3.aaToMat([1 0 0],pi/2);
		E2(1:3,4) = [-l/2-radius 0 0]';
		model.constraints{end+1} = apbd.ConJointRev({model.bodies{end},model.bodies{end-1}},{E2,E1});
		model.constraints{end}.radius = radius*2;
		model.constraints{end}.height = 1.5;
        model.bodies{end}.setInitVelocity([0 -0 0 0 0 0]');

		n = 4;
        sides = [w w w];
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
			E = eye(4);
			x = 0.05*i-0.4;
			y = 0;
			z = (i-0.5)*w*0.99;
			E(1:3,4) = [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end

    	case 2
		model.name = 'Rigid Collisions: Pyramid';
		model.plotH = false;
		model.tEnd = 1;
		model.h = 5e-3;
		model.substeps = 1;
		model.iters = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.0;

		model.ground.size = 10;
		model.axis = 5*[-1 1 -1 1 0 1];
		model.drawHz = 10000;

		model.view = [0 0];

		layers = 3;
        for i = 1: layers
		    for j = 1 : i
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
			    E = eye(4);
			    x = -0.5 *w *(i-1) + j* w;
			    y = 0;
			    z = (layers - i + 0.5)*w*0.99;
			    E(1:3,4) = [x y z]';
			    model.bodies{end}.setInitTransform(E);
                if i == 2 && j == 1
                    model.bodies{end}.setInitVelocity([0 0 0 10 0 0]');
                end
            end
        end
end

end
