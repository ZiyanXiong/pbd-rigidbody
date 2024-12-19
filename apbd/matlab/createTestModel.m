function model = createTestModel(modelID, h, substeps)

model = apbd.Model();

switch(modelID)
    	case 0
		model.name = 'Stacking: 10 rigid bodies with offset';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.0;

		model.ground.size = 20;
		model.axis = 2*[-1 1 -1 1 0 1];
		model.drawHz = 1000;

		model.view = [0 0];

		n = 2;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  0.1 * i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
    	case 1
		model.name = 'Stacking: 10 rigid bodies with offset';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
    	case 2
		model.name = 'Stacking: 2 rigid bodies without fricition';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.0;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 2;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.5*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
    	case 3
		model.name = 'Dynamic Fricition: sliding distance test';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.9;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 1;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 100 0 0]');
            end
        end
	case 4
		model.name = 'Static Friciton: 2 rigid bodies on a slope';
		model.plotH = false;
		model.tEnd = 1.0;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 20;
		density = 1.0;
		l = 1;
		w = 1;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		angle = 20*pi/180;
		R = se3.aaToMat([0 1 0],angle);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;
        
        n = 2;
        for i = 1 : n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = 1.01*(sin(angle)/cos(angle));
		    E = eye(4);
		    E(1:3,1:3) = R;
			x = 0.00*i;
			y = 0;
			z = (i-0.5)*w;
			E(1:3,4) = R* [x y z]';
		    model.bodies{end}.setInitTransform(E);
            if i == 2
		        model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
	case 5
		model.name = 'Static Friciton: 10 rigid bodies on a slope';
		model.plotH = false;
		model.tEnd = 1.0;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 20;
		density = 1.0;
		l = 1;
		w = 1;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		angle = 20*pi/180;
		R = se3.aaToMat([0 1 0],angle);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;
        
        n = 10;
        for i = 1 : n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = 1.01*(sin(angle)/cos(angle));
		    E = eye(4);
		    E(1:3,1:3) = R;
			x = 0.00*i;
			y = 0;
			z = (i-0.5)*w;
			E(1:3,4) = R* [x y z]';
		    model.bodies{end}.setInitTransform(E);
            if i == 2
		        model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
    	case 6
		model.name = 'Stacking: 10 rigid bodies falling one by one';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i;
			y = 0;
			z = (i-0.5 + i * 0.1)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
	case 7
		model.name = 'Stacking: 10 rigidbodies and push the second one';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 100 0 0]');
            end
        end
        case 8
		model.name = 'Stacking: 10 rigidbodies with large offset';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.4*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        case 9
		model.name = 'Stacking: 10 rigidbodies tetris';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w + 1;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        case 10
		model.name = 'Stacking: the wall';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
        l = 2;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		for i = 1 : 6
            for j = 1 : 3
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], 0.0);
			    E = eye(4);
			    x = mod(i+1,2)*0.45*l + (j-0.5) * l - l*2.5 + j * 0.0;
			    y = 0;
			    z = (i-0.5 + i *0.0)*w;
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
                if i == 2
                    model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
                end
            end
        end
        case 11
		model.name = 'Stacking: heavy head';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
            model.bodies{end}.density = 2^((i-1)) * 1;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        case 12
		model.name = 'Stacking: heavy head';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 35*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 5;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(2^(i-1) * sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.0*i;
			y = 0;
			z = 0.5*w*2^(i-1) + 2^(i-1) - 1;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        case 13
		model.name = 'Stacking: seasaw';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 10;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 6;
        for j = 1 : 2
		    for i = 1 : n
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], 0.0);
			    E = eye(4);
                if (j == 1)
			        x = 1.5;
                else
                    x = -1.5;
                end
			    y = 0;
			    z = (i-0.5)*w + 2;
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
                if i == 2 && j == 1 
                    break;
                end
            end
        end
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
        E = eye(4);
        E(3,4) = 0.5;
        model.bodies{end}.setInitTransform(E);
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([4 1 1]'),density);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
        E = eye(4);
        E(3,4) = 1.5;
        model.bodies{end}.setInitTransform(E);
        case 14
		model.name = 'Stacking: seasaw';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 10;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 6;
        for j = 1 : 2
		    for i = 1 : n
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], 0.0);
			    E = eye(4);
                if (j == 1)
			        x = 1.5;
                else
                    x = -1.5;
                end
			    y = 0;
			    z = (i-0.5)*w + 2;
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
            end
        end
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
        E = eye(4);
        E(3,4) = 0.5;
        model.bodies{end}.setInitTransform(E);
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([4 1 1]'),density);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
        E = eye(4);
        E(3,4) = 1.5;
        model.bodies{end}.setInitTransform(E);
        case 16
		model.name = 'Stacking: 10 rigidbodies tetris';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 2;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 4;
        xs = 2*[-1 -1  1  1];
        ys = 2*[-1  1 -1  1];
		for i = 1 : n
            for j = 1 : 1
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([1 1 4]),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], 0.0);
			    E = eye(4);
			    x = xs(i);
			    y = ys(i);
			    z = 2+ 5* (j-1);
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
            end
        end
        
		for i = 1 : 1
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([6 6 1]),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
		    E = eye(4);
		    x = 0;
		    y = 0;
		    z = -0.5 + 5*i;
            E(1:3,1:3) = R;
		    E(1:3,4) = R * [x y z]';
		    model.bodies{end}.setInitTransform(E);
        end
    	case 17
		model.name = 'Stacking: 10 mesh rigid bodies with offset';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 2;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1]*w;
		model.drawHz = 60;

		model.view = [0 0];
        meshShape = apbd.ShapeMeshObj('./ShapeFiles/cube.obj');
        meshShape.computeInertia(density);
		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.04*i*w;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
    	case 18
		model.name = 'Stacking: Mesh';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = substeps;
        %model.itersSP = 30;
		density = 1.0;
		w = 1;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 2*w*[-1 1 -1 1 0 1];
		model.drawHz = 1000;

		model.view = [0 0];

		angle = -90*pi/180;
        meshShape = apbd.ShapeMeshObj('./ShapeFiles/bunny.obj');
        meshShape.computeInertia(density);
		n = 2;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            %R = se3.aaToMat([0 0 1], 0.0);
            R = se3.aaToMat([0 0 1],angle);
			E = eye(4);
			x = 0 * w * (i-1);
			y = 0.75 * w * (i-2);
			%z = (i-0.5 + i *0.0)*w;
            z = -0.15 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E*meshShape.E_oi);
            if i == 2
                %model.bodies{end}.setInitVelocity([0 0 0 1 0 0]');
            end
        end
        
		model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        %R = se3.aaToMat([0 0 1], 0.0);
        R = se3.aaToMat([0 0 1],angle);
		E = eye(4);
		x = -1;
		y = -0.75;
		%z = (i-0.5 + i *0.0)*w;
        z = -0.15 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
        model.bodies{end}.setInitTransform(E*meshShape.E_oi);

        
		model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        %R = se3.aaToMat([0 0 1], 0.0);
        R = se3.aaToMat([0 0 1],0);
		E = eye(4);
		x = -0.5;
		y = 0.0 * w * (i-2);
		%z = (i-0.5 + i *0.0)*w;
        z = 1.25 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
        model.bodies{end}.setInitTransform(E);
    	case 19
		model.name = 'Stacking: Cube Mesh';
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = substeps;
        %model.itersSP = 30;
		density = 1.0;
		w = 2;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 2*w*[-1 1 -1 1 0 1];
		model.drawHz = 1000;

		model.view = [0 0];

		angle = -90*pi/180;
        meshShape = apbd.ShapeMeshObj('./ShapeFiles/cube.obj');
        meshShape.computeInertia(density);
		n = 2;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            %R = se3.aaToMat([0 0 1], 0.0);
            R = se3.aaToMat([0 0 1],angle);
			E = eye(4);
			x = 0 * w * (i-1);
			y = 1.25 * w * (i-2);
			%z = (i-0.5 + i *0.0)*w;
            z = -0.0 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E*meshShape.E_oi);
            if i == 2
                %model.bodies{end}.setInitVelocity([0 0 0 1 0 0]');
            end
        end
        
		model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        %R = se3.aaToMat([0 0 1], 0.0);
        R = se3.aaToMat([0 0 1],angle);
		E = eye(4);
		x = -2.1;
		y = -1.25;
		%z = (i-0.5 + i *0.0)*w;
        z = -0.0 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
        model.bodies{end}.setInitTransform(E*meshShape.E_oi);
        
		model.bodies{end+1} = apbd.BodyRigid(meshShape, density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        %R = se3.aaToMat([0 0 1], 0.0);
        R = se3.aaToMat([0 0 1],0);
		E = eye(4);
		x = -1.25;
		y = 1;
		%z = (i-0.5 + i *0.0)*w;
        z = 1.51 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
        model.bodies{end}.setInitTransform(E);
        case 20
		model.name = 'Stacking: Catenary';
		model.plotH = false;
		model.tEnd = 5;
		model.h = h;
		model.substeps = substeps;
		model.iters = substeps;
        %model.itersSP = 3;
		density = 1.0;
		w = 3;
		
		model.grav = [0 0 -981]';
		model.ground.E = eye(4);
		mu = 0.9;

		model.ground.size = 10;
		model.axis = 20*[-1 1 -1 1 0 2.5];
		model.drawHz = 60;

		model.view = [0 0];

        a = 6.4;
		n = 9;

        xs = [-15 -12 -9 -5 0 5 9 12 15];
        halfAngles = [7.5 2.5 5 20 20 20 5 2.5 7.5] ./ 180 * pi;
        %halfAngles = ones(9,1) * (10/180) * pi;
        halfDistances = [1.1 0.7 0.5 0.5 0.3 0.5 0.5 0.7 1.1]* w;
        centers = zeros(3,n);
        lengths = zeros(n,1);
        heights = [w*1.6 w*1.5 w*1.2 w*0.8 w*0.6 w*0.8 w*1.2 w*1.5 w*1.6];
        densities = [0.1 0.2 0.4 0.6 1 0.6 0.4 0.2 0.1];
        %densities = ones(9,1);
        theta = 0;
        for i = 1:n
            centers(1,i) = xs(i);
            centers(2,i) = 0;
            centers(3,i) = 35-(a*cosh(centers(1,i)/a)-a);      
            if(i==1)
                lengths(i) = 2*(centers(3,i) - halfDistances(i) * cos(halfAngles(i)));
            else
                lengths(i) = 2*(abs((centers(:,i) - centers(:,i-1))'*[cos(pi/2 - theta) 0 sin(pi/2 - theta)]') - lengths(i-1)/2 ...
                    - halfDistances(i-1) * cos(halfAngles(i-1)) - halfDistances(i)* cos(halfAngles(i)));
            end
            theta = theta + 2*halfAngles(i);  
        end
        lengths(5) = lengths(5);
        theta = pi/2;
		for i = 1 : n
            sides = [lengths(i) w heights(i)];
            model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeTwoCuboid(sides, sides, halfDistances(i), halfAngles(i)),densities(i));
            %model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
            
		    R = se3.aaToMat([0 1 0], theta + halfAngles(i));
			E = eye(4);
            E(1:3,1:3) = R;
			E(1:3,4) = centers(:,i);
			model.bodies{end}.setInitTransform(E);
            theta = theta + 2*halfAngles(i);
        end 
    	case 21
		model.name = 'Stacking : Arch';
		model.plotH = false;
		model.tEnd = 5;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 3;
		density = 1.0;
		w = 3;
		sides = [w w w];
		model.grav = [0 0 -981]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 10;
		model.axis = 20*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
        halfAngle = 0.5 * pi / n;
        halfDistance = 0.4 * w;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeTwoCuboid(sides, sides, halfDistance, halfAngle),density);
            %model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
            theta = (i*2-1)*halfAngle;
            r = (0.5*w + cos(halfAngle) * halfDistance) / sin(halfAngle);
    		R = se3.aaToMat([0 1 0], pi/2 + theta);
			E = eye(4);
			x = -r * cos(theta);
			y = 0;
			z = r*sin(theta);
            E(1:3,1:3) = R;
			E(1:3,4) = [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                %model.bodies{end}.setInitVelocity([0 0 0 0 0 0]', model.h);
            end
        end
end

end
