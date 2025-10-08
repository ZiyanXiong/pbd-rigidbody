function model = JointModel(modelID, h, substeps, solverType)

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
		model.axis = 10*w*[-1 1 -1 1 0 1];
		model.drawHz = 10;

		model.view = [0 0];
        model.solverType = solverType;
        groundCollisionList = [];
        bodyCollisionList = [];
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
            if(i == 2)
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
            groundCollisionList(end+1) = i;
            if(i~=n)
                bodyCollisionList(end+1,:) = [i i+1];
            end
        end

        %bodyCollisionList = [];
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
        
    case 1
		model.name = 'Joint:2 Hinge Joint';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 5;
		model.h = h;
		model.steps = ceil(model.tEnd/model.h);
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1;
		w = 4;
		sides = [w 2.5*w w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = 2.5 * w *[-5 5 -1 4 0 5];
		model.drawHz = 10;

		model.view = [90 0];

	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),Inf);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([1 0 0], 0);
	    E = eye(4);
	    x = 0;
	    y = 1.25*w;
	    z = 0;
        E(1:3,1:3) = R;
	    E(1:3,4) = R * [x y z]' + [0 0 10.5*w]';
	    model.bodies{end}.setInitTransform(E);

        n = 2;
        for i = 1:n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            if(i~=n)
                R = se3.aaToMat([1 0 0], 0);
                E = eye(4);
                x = 0;
                y = (2.5*i+1.25)*w;
                z = 0;
                E(1:3,1:3) = R;
                E(1:3,4) = R*[x y z]' + [0 0 10.5*w]';
            else
                R = se3.aaToMat([1 0 0], 0);
                E = eye(4);
                x = 0;
                y = (2.5*i+1.25)*w;
                z = 0;
                E(1:3,1:3) = R;
                E(1:3,4) = R*[x y z]' + [0 0 10.5*w]';
            end
		    model.bodies{end}.setInitTransform(E);
            model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
        end
        
        m = 10;
        for i = 1:m
	        model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid([w w w]),density);
	        model.bodies{end}.collide = true;
	        model.bodies{end}.mu = mu;
	        %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            R = se3.aaToMat([1 0 0], 0);
	        E = eye(4);
	        x = 0;
	        y = 6*w;
	        z = (i-0.5)*w;
            E(1:3,1:3) = R;
	        E(1:3,4) = R * [x y z]';
	        model.bodies{end}.setInitTransform(E);
        end
        
        for i = 1:n
            model.joints{end+1} = JointHinge(model.bodies{i}, model.bodies{i+1}, false, [0 2.5*w*i 10.5*w]' ,[1 0 0]', -000000*(6-i)*ones(model.steps,1));
        end
        groundCollisionList = n+2:n+m+1;
        bodyCollisionList = [];
        for i = 1:m-1
            bodyCollisionList(end+1,:) = [i+n+1 i+n+2];
        end
        for i = 1:m
            bodyCollisionList(end+1,:) = [n+1 n+1+i];
        end

        %bodyCollisionList = [];
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
    	case 2
		model.name = 'Stacking : Arch';
        model.modelID = modelID;
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
		model.axis = 30*[-1 1 -1 1 0 1];
		model.drawHz = 60;

		model.view = [0 0];
        model.solverType = solverType;
        groundCollisionList = [];
        bodyCollisionList = [];

		n = 3;
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
            groundCollisionList(end+1) = i;
            if(i~=n)
                bodyCollisionList(end+1,:) = [i i+1]';
            end
        end
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
        case 3
		model.name = 'Joint:Torque on joints';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 2;
		model.h = h;
		model.steps = ceil(model.tEnd/model.h);
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1;
		w = 6;
		sides = [w 1 1];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = w *[-2 2 -2 2 0 4];
		model.drawHz = 10;

		model.view = [0 0];
        
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),inf);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([1 0 0], 0);
	    E = eye(4);
	    x = 0;
	    y = 0;
	    z = 2.5*w;
        E(1:3,1:3) = R;
	    E(1:3,4) = [x y z]';
	    model.bodies{end}.setInitTransform(E);
        
        n = 2;
        for i = 1:n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
	        %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            if(i==1)
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*(i-0.5);
                y = 0;
                z = 2.5*w;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]'  + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            else
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*(i-0.5);
                y = 0;
                z = 2.5*w ;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]' + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 100 0 0 0]');
            end
        end
        %{
        torques = -100 * ones(model.steps,1);
        load("sin_torques_2.mat");
        torques1 = torques1 / h;
        torques2 = torques2 / h;
        timesteps = 1:model.steps;
        targets1 = 0.5*sin(timesteps * 2*pi/ model.steps);
        targets2 = 0.5*sin(timesteps * 4*pi/ model.steps);
        %ts(1,1) = -100000;
        %}
        targets1 = pi/4 * ones(model.steps,1);
        targetsW1 = 0 * ones(model.steps,1);
        targets2 = pi/4 * ones(model.steps,1);
        targetsW2 = 0 * ones(model.steps,1);
        for i = 1:n
            if(i==1)
                %model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', torques1);
                model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', targets1, targetsW1, 1000, 10);
                %model.joints{end+1} = JointFix(model.bodies{i}, model.bodies{i+1}, false, [0 0 11.5*w]', [0 2.5*w*(i-1)+2*w 10.5*w]', norm([0 0 11.5*w]'- [0 2.5*w*(i-1)+2*w 10.5*w]'));
            else
                %model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', torques2);
                model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', targets2, targetsW2, 1000, 10);
            end
        end
        groundCollisionList = 1:2;
        bodyCollisionList = [];
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
    case 4
		model.name = 'Joint:Free joints';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 2;
		model.h = h;
		model.steps = ceil(model.tEnd/model.h);
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1;
		w = 6;
		sides = [w 1 1];
		model.grav = [0 0 0]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = w *[-2 2 -2 2 0 4];
		model.drawHz = 10;

		model.view = [0 0];
        
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),inf);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([1 0 0], 0);
	    E = eye(4);
	    x = 0;
	    y = 0;
	    z = 2.5*w;
        E(1:3,1:3) = R;
	    E(1:3,4) = [x y z]';
	    model.bodies{end}.setInitTransform(E);
        
        n = 2;
        for i = 1:n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
	        %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            if(i==1)
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*(i-0.5);
                y = 0;
                z = 2.5*w;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]'  + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            else
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*(i-0.5);
                y = 0;
                z = 2.5*w ;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]' + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 10 0 0 0]');
            end
        end
        torques = zeros(model.steps,1);
        for i = 1:n
            if(i==1)
                model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', torques);
                %model.joints{end+1} = JointSpring(model.bodies{i}, model.bodies{i+1}, false, [0 0 0]' ,[0 0 0]', 4, 50000);
                %model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', targets1, targetsW1, 1000, 10);
                %model.joints{end+1} = JointFix(model.bodies{i}, model.bodies{i+1}, false, [0 0 11.5*w]', [0 2.5*w*(i-1)+2*w 10.5*w]', norm([0 0 11.5*w]'- [0 2.5*w*(i-1)+2*w 10.5*w]'));
            else
                model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', torques);
                model.joints{end+1} = JointSpring(model.bodies{i}, model.bodies{i+1}, false, [0 0 0]' ,[0 0 0]', 5, 5000);
                %model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', targets2, targetsW2, 1000, 10);
            end
        end

        groundCollisionList = 1:2;
        bodyCollisionList = [];
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
    case 5
		model.name = 'Joint:Spring';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 2;
		model.h = h;
		model.steps = ceil(model.tEnd/model.h);
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1;
		w = 4;
		sides = [w w w];
		model.grav = [0 0 0]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = w *[-2 2 -2 2 0 4];
		model.drawHz = 10;

		model.view = [0 0];

        n = 2;
        for i = 1:n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
	        %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            if(i==1)
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*1.5;
                y = 0;
                z = 0.5*w;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]'  + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            else
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = -w*1.5;
                y = 0;
                z = 0.5*w ;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]' + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        torques = zeros(model.steps,1);
        for i = 1:n
            if(i==1)
                %model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', torques);
                %model.joints{end+1} = JointSpring(model.bodies{i}, model.bodies{i+1}, false, [-2 0 0]' , [2 0 0]', 5, 5000);
                %model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', targets1, targetsW1, 1000, 10);
                %model.joints{end+1} = JointFix(model.bodies{i}, model.bodies{i+1}, false, [0 0 11.5*w]', [0 2.5*w*(i-1)+2*w 10.5*w]', norm([0 0 11.5*w]'- [0 2.5*w*(i-1)+2*w 10.5*w]'));
                muscleBodies = {model.bodies{i}, model.bodies{i+1}};
                points = [[0 0 0]',[-2 0 0]', [2 0 0]', [0 0 0]'];
                model.muscles{end+1} = MuscleSpring(muscleBodies, false, points, 9, 5000);
            else
                %model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', torques);
                %model.joints{end+1} = JointSpring(model.bodies{i}, model.bodies{i+1}, false, [0 0 0]' ,[0 0 0]', 5, 1000);
                %model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', targets2, targetsW2, 1000, 10);
            end
        end

        groundCollisionList = [];
        bodyCollisionList = [];
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
    case 6
		model.name = 'Muscle: Via points';
        model.modelID = modelID;
		model.plotH = false;
		model.tEnd = 2;
		model.h = h;
		model.steps = ceil(model.tEnd/model.h);
		model.substeps = substeps;
		model.iters = 1;
        model.solverType = solverType;

        %model.itersSP = 30;
		density = 1;
		w = 6;
		sides = [w 1 1];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.5;

		model.ground.size = 20;
		model.axis = w *[-2 2 -2 2 0 4];
		model.drawHz = 10;

		model.view = [0 0];
        
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),inf);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([1 0 0], 0);
	    E = eye(4);
	    x = -1.5*w;
	    y = 0;
	    z = 3*w;
        E(1:3,1:3) = R;
	    E(1:3,4) = [x y z]';
	    model.bodies{end}.setInitTransform(E);
        
        n = 3;
        for i = 1:n
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
	        %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            if(i==1)
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*(i-2);
                y = 0;
                z = 3*w;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]'  + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            else
                R = se3.aaToMat([0 1 0], 0);
                E = eye(4);
                x = w*(i-2);
                y = 0;
                z = 3*w ;
                E(1:3,1:3) = R;
                E(1:3,4) = R *[0.5*w 0 0]' + [x y z]';
		        model.bodies{end}.setInitTransform(E);
                model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
            end
        end
        torques = zeros(model.steps,1);
        for i = 1:n
            if(i==1)
                model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', torques);
                model.joints{end}.setLimits(pi/4, -pi/4);
                %model.joints{end+1} = JointSpring(model.bodies{i}, model.bodies{i+1}, false, [0 0 0]' ,[0 0 0]', 4, 50000);
                %model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', targets1, targetsW1, 1000, 10);
                %model.joints{end+1} = JointFix(model.bodies{i}, model.bodies{i+1}, false, [0 0 11.5*w]', [0 2.5*w*(i-1)+2*w 10.5*w]', norm([0 0 11.5*w]'- [0 2.5*w*(i-1)+2*w 10.5*w]'));
                %muscleBodies = {model.bodies{2}, model.bodies{3}, model.bodies{4}};
                %points = [[0 0 0]',[0 0 -0.5]', [0 0 -0.5]', [2 0 -0.5]', [0 0 0]'];
                %model.muscles{end+1} = MuscleSpring(muscleBodies, false, points, 12, 5000);
            else
                model.joints{end+1} = JointHinge2(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 1 0]', torques);
                model.joints{end}.setLimits(pi/4, -1/4*pi);
                %model.joints{end+1} = JointSpring(model.bodies{i}, model.bodies{i+1}, false, [0 0 0]' ,[0 0 0]', 5, 5000);
                %model.joints{end+1} = JointHinge2Actuator(model.bodies{i}, model.bodies{i+1}, false, [0.5*w 0 0]' ,[0 0 1]', targets2, targetsW2, 1000, 10);
            end
        end

        groundCollisionList = 1:2;
        bodyCollisionList = [];
        model.collider = apbd.Collider(model,groundCollisionList, bodyCollisionList);

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
            fid = fopen(fullfile(model.resultFolder, 'Body_States_GPQP.txt'), 'w');
            fclose(fid);
        end
end
