function model = SIG25Model(modelID, h, substeps, solverType)

model = apbd.Model();

switch(modelID)
    case 0
		model.name = 'Scene:Test';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 10*[-1 1 -1 1 0 1];
		model.drawHz = 10;

		model.view = [0 0];
        model.solverType = solverType;

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  0.0 * i;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
    case 1
		model.name = 'Stacking: 20 rigid bodies vertically';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 3;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 1];
		model.drawHz = 120;

		model.view = [0 0];
        model.solverType = solverType;

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x =  0.0;
		y = 0;
		z = 0.5*w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);
        angle = 0;
		n = 3;
		for i = 2 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 1 0], angle);
			E = eye(4);
			x =  0.2 ;
			y = 0;
			z = (i-1.5)*w+0.1;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]' + [-0.5*w*sin(angle) 0 w+0.5*w*sin(angle)]';
			model.bodies{end}.setInitTransform(E);
            model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
    case 2
		model.name = 'Stacking: 10 rigid bodies vertically with noise';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 1];
		model.drawHz = 120;

		model.view = [0 0];
        model.solverType = solverType;

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], rand * pi* 0.05);
			E = eye(4);
			x =  (rand - 0.5)*0.05;
			y = (rand - 0.5)*0.05;
			z = (i-0.5 + i *0.0)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
    case 3
		model.name = 'Static Friciton: 10 rigid bodies on a slope';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1.0;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 20;
		density = 1.0;
		l = 4;
		w = 4;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		angle = 15*pi/180;
		R = se3.aaToMat([0 1 0],angle);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 20;
		model.axis = 60*[-1 1 -1 1 0 2];
		model.drawHz = 60;
        model.solverType = solverType;
        
        n = 10;
        for i = 1 : n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = 1.5*(sin(angle)/cos(angle));
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
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
	case 4
		model.name = 'Static Friciton: 10 rigid bodies on a slope (verical)';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1.0;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 20;
		density = 1.0;
		l = 4;
		w = 4;
		sides = [l w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);

		angle = 15*pi/180;
		R = se3.aaToMat([0 1 0],angle);
		model.ground.E(1:3,1:3) = R;

		model.ground.size = 20;
		model.axis = 60*[-1 1 -1 1 0 2];
		model.drawHz = 60;
        model.solverType = solverType;
        
        n = 10;
        for i = 1 : n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = 0.5;
		    E = eye(4);
		    E(1:3,1:3) = R;
			x = -0.25*w*i;
			y = 0;
			z = (i-0.5)*w;
			E(1:3,4) = R* [x y z]';
		    model.bodies{end}.setInitTransform(E);
            if i == 2
		        model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
	case 5
	    model.name = 'Stacking: Inverted tower';
        model.modelID = modelID;
	    model.plotH = false;
	    model.tEnd = 1;
	    model.h = h;
	    model.substeps = substeps;
	    model.iters = 1;
        %model.itersSP = 30;
	    density = 1.0;
	    w = 4;
	    sides = [w w w];
	    model.grav = [0 0 -980]';
	    model.ground.E = eye(4);
	    mu = 0.5;
    
	    model.ground.size = 20;
	    model.axis = 80*[-1 1 -1 1 0 2];
	    model.drawHz = 120;
    
	    model.view = [0 0];
        model.solverType = solverType;
    
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
			z = 0.5*w*2^(i-1) + (2^(i-1) - 1)*w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 2
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
    case 6
		model.name = 'Stacking: Door';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 2];
		model.drawHz = 120;

		model.view = [0 0];
        model.solverType = solverType;

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  -5*w;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end

		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  5*w;
			y = 0;
			z = (i-0.5 + i *0.0)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([11*w w w]),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x =  0;
		y = 0;
		z = (n+0.5)*w;
        %z = 0.5 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);

        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
	case 7
		model.name = 'Stacking: Balanced goal post';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 2];
		model.drawHz = 120;

		model.view = [0 0];
        model.solverType = solverType;


		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x =  0;
		y = 0;
		z = 0.5*w;
        %z = 0.5 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);
        
		n = 9;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  -5*w;
			y = 0;
			z = (i-0.5 + 2)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end

		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  5*w;
			y = 0;
			z = (i-0.5 + 2)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end
        
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([11*w w w]),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x =  0;
		y = 0;
		z = 1.5*w;
        %z = 0.5 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);

        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 8
		model.name = 'Stacking: Unbalanced goal post';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 2];
		model.drawHz = 120;

		model.view = [0 0];
        model.solverType = solverType;


		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x =  0;
		y = 0;
		z = 0.5*w;
        %z = 0.5 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);
        
		n = 8;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  -5*w;
			y = 0;
			z = (i-0.5 + 2)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end

		for i = 1 : n + 1
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x =  5*w;
			y = 0;
			z = (i-0.5 + 2)*w;
            %z = 0.5 * w;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
        end
        
		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([11*w w w]),density);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x =  0;
		y = 0;
		z = 1.5*w;
        %z = 0.5 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density*2.7);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x = -5*w;
		y = 0;
		z = (n-0.5 + 4.0)*w;
        %z = 0.5 * w;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);

        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = true;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 9
		model.name = 'Stacking: EarthQuake';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
        model.solverType = solverType;

		density = 1.0;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 60*[-1 1 -1 1 0 2];
		model.drawHz = 10;

		model.view = [0 0];

		model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([80 80 0.2]), Inf);
		model.bodies{end}.collide = true;
		model.bodies{end}.mu = mu;
		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0.0);
		E = eye(4);
		x = 0.0;
		y = 0;
		z = 1;
        E(1:3,1:3) = R;
		E(1:3,4) = R * [x y z]';
		model.bodies{end}.setInitTransform(E);

		n = 10;
		for i = 1 : n
			model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			model.bodies{end}.collide = true;
			model.bodies{end}.mu = mu;
    		%R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([0 0 1], 0.0);
			E = eye(4);
			x = 0.0*i;
			y = 0;
			z = (i-0.5 + i *0.0)*w + 1.1;
            E(1:3,1:3) = R;
			E(1:3,4) = R * [x y z]';
			model.bodies{end}.setInitTransform(E);
            if i == 1
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = false;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 10
		model.name = 'Stacking: Jenga Add';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        %model.itersSP = 30;
        model.solverType = solverType;

		density = 0.6;
		w = 4;
		sides = [2*w 6*w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.2;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 1];
		model.drawHz = 10;

		model.view = [45 45];

        layers = 10;
        for l = 1:layers
		    for i = 1 : 3
                if((l==4||l==9)&&(i==3))
                    continue;
                end
                if((l==10)&&(i==2||i==1))
                    continue;
                end
                if((l<4)&&(i==3||i==1))
                    continue;
                end
                if((l<8&&l>5)&&(i==1))
                    continue;
                end
                
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], pi/2 * mod(l+1,2));
			    E = eye(4);
			    x = -4.03*w + 2.01*w*i;
			    y = 0;
			    z =-0.5*w + w*l;
                if(l==10)
                    x = x - 3.5 * w;
                    y = y + 3.5*w;
                    z = z + 1*w;
                end
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
            end
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            model.useContactCaching = false;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 11
		model.name = 'Stacking: Jenga Remove';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 0.6;
		w = 4;
		sides = [2*w 6*w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.2;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 1];
		model.drawHz = 10;

		model.view = [0 0];

        layers = 10;
        for l = 1:layers
		    for i = 1 : 3
                if((l==4||l==9)&&(i==3))
                    continue;
                end
                if((l==10)&&(i==2||i==1))
                    continue;
                end
                if((l<4)&&(i==3||i==1))
                    continue;
                end
                if((l<8&&l>5)&&(i==1))
                    continue;
                end
                
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], pi/2 * mod(l+1,2));
			    E = eye(4);
			    x = -4.03*w + 2.01*w*i;
			    y = 0;
			    z =-0.5*w + w*l;
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
            end
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = false;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 12
		model.name = 'Stacking: Bowls';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		%sides = [2*w 6*w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.6;

		model.ground.size = 20;
		model.axis = 30*[-1 1 -1 1 0 1];
		model.drawHz = 10;

		model.view = [0 0];
        fileNames = {'./ShapeFiles/bowl/bowl.obj', ...
            './ShapeFiles/bowl/bowl_part2.obj', ...
            './ShapeFiles/bowl/bowl_part3.obj', ...
            './ShapeFiles/bowl/bowl_part4.obj', ...
            './ShapeFiles/bowl/bowl_part5.obj', ...
            './ShapeFiles/bowl/bowl_part6.obj'};

        n = 10;
        for i = 1:n
	        model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeNonConvexObj(fileNames),density);
	        model.bodies{end}.collide = true;
	        model.bodies{end}.mu = mu;
	        %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
	        E = eye(4);
	        x = 0;
	        y = 0;
	        z =(1.0*i-0.5)*w*0.49 + 0.3;
            if(mod(i,2)==0)
                R = se3.aaToMat([0 1 0], 3 * pi/180);
                x = -0.03 * w;
            else
                x = 0.03 * w;
                R = se3.aaToMat([0 1 0], -3 * pi/180);
            end
            if(i == 1)
                x = 0;
                R = eye(3);
            end
            E(1:3,1:3) = R;
	        E(1:3,4) = [x y z]';
	        model.bodies{end}.setInitTransform(E);
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = false;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 13
		model.name = 'Stacking: Circular';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1.0;
		w = 4;
		sides = [2*w 3.5*w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 20*[-1 1 -1 1 0 2];
		model.drawHz = 10;

		model.view = [0 0];

        layers = 10;
        for l = 1:layers
		    for i = 1 : 5                
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], i*2*pi/5 + pi/6 * l);
			    E = eye(4);
			    x = -4*w;
			    y = 0;
			    z =-0.5*w + w*l;
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
            end
        end
        
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([0 0 1], 0);
	    E = eye(4);
	    x = -0.0 *w;
	    y = 14*w;
	    z =-0.5*w + w*7;
        E(1:3,1:3) = R;
	    E(1:3,4) = R * [x y z]';
	    model.bodies{end}.setInitTransform(E);
        model.bodies{end}.setInitVelocity([0 0 0 0 -600 150]');
        
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = false;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
        case 14
		model.name = 'Stacking: Jenga Test';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 1;
		model.h = h;
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 0.6;
		w = 4;
		sides = [2*w 6*w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.2;

		model.ground.size = 20;
		model.axis = 40*[-1 1 -1 1 0 1];
		model.drawHz = 10;

		model.view = [0 0];

        layers = 3;
        for l = 1:layers
		    for i = 1 : 2
			    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
			    model.bodies{end}.collide = true;
			    model.bodies{end}.mu = mu;
    		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
                R = se3.aaToMat([0 0 1], pi/2 * mod(l+1,2));
			    E = eye(4);
			    x = -4.03*w + 2.01*w*i;
			    y = 0;
			    z =-0.5*w + w*l;
                E(1:3,1:3) = R;
			    E(1:3,4) = R * [x y z]';
			    model.bodies{end}.setInitTransform(E);
            end
        end
        model.resultFolder = sprintf("Results\\Scene\\%d\\",model.modelID);
        if ~exist(model.resultFolder, 'dir')
           mkdir(model.resultFolder)
        end
        if ~exist(strcat(model.resultFolder,"residual_per_iteration\\"), 'dir')
           mkdir(strcat(model.resultFolder,"residual_per_iteration\\"))
        end
        if(model.solverType == 1)
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        elseif(model.solverType == 2)
            fid = fopen(fullfile(model.resultFolder, 'Body_States_2PSP.txt'), 'w');
            fclose(fid);
        else
            %model.useContactCaching = false;
            fid = fopen(fullfile(model.resultFolder, sprintf('Body_States_TGS_%d.txt',model.substeps)), 'w');
            fclose(fid);
        end
end

end
