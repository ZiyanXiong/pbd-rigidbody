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
			x = i - 1 + 0.1 * i;
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
		w = 1;
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
		model.name = 'Two Cuboid Rigid Collisions : arch';
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
		mu = 0.9;

		model.ground.size = 10;
		model.axis = 20*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];

		n = 10;
        halfAngle = (9/180) * pi;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeTwoCuboid(sides, sides, 0.4*3, halfAngle),density);
            %model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
            theta = (i*2-1)*halfAngle;
            r = 0.89 / sin(halfAngle) * 3;
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
        %{
		for i = 1 : 3
			%model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeTwoCuboid(sides, sides, 0.4, halfAngle),density);
            model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
			E = eye(4);
			x = 0.05*i;
			y = 0;
			z = (i-0.5)*w+6.3;
			E(1:3,4) = [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        %}
        %model.constraintList{end+1} = apbd.ConFixLocalPoint(model.bodies{end}, [-0.5,-0.5,0-0.5]');
        %model.constraintList{end+1} = apbd.ConFixLocalPoint(model.bodies{end}, [0.5,0.5,0-0.5]');
        %model.constraintList{end+1} = apbd.ConFixLocalPoint(model.bodies{end}, [-0.5,0.5,0-0.5]');
        %model.constraintList{end+1} = apbd.ConFixLocalPoint(model.bodies{end}, [0.5,-0.5,0-0.5]');
end

end
