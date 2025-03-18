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
		model.axis = 10*[-1 1 -1 1 0 1];
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
                model.bodies{end}.setInitVelocity([0 0 0 100 0 0]');
            end
            groundCollisionList(end+1) = i;
            if(i~=n)
                bodyCollisionList(end+1,:) = [i i+1]';
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
		sides = [2*w 6*w 2*w];
		model.grav = [0 0 -980]';
		model.ground.E = eye(4);
		mu = 0.2;

		model.ground.size = 20;
		model.axis = 2.5 * w *[-5 5 -5 5 0 10];
		model.drawHz = 10;

		model.view = [90 0];

	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),Inf);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([1 0 0], pi/2);
	    E = eye(4);
	    x = 0;
	    y = 3*w;
	    z = 0;
        E(1:3,1:3) = R;
	    E(1:3,4) = R * [x y z]' + [0 0 10*w]';
	    model.bodies{end}.setInitTransform(E);

        n = 3;
        for i = 1:n-1
		    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
		    model.bodies{end}.collide = true;
		    model.bodies{end}.mu = mu;
		    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
            if(i<5)
                R = se3.aaToMat([1 0 0], 0);
	            E = eye(4);
	            x = 0;
	            y = -3*w + i*6*w;
	            z = 0;
                E(1:3,1:3) = R;
	            E(1:3,4) = R*[x y z]' + [0 0 10*w]';
            else
                R = se3.aaToMat([1 0 0], pi/2);
	            E = eye(4);
	            x = 0;
	            y = 3*w;
	            z = 0;
                E(1:3,1:3) = R;
	            E(1:3,4) = R * [x y z]' + [0 0 i*6*w]';
            end
		    model.bodies{end}.setInitTransform(E);
            model.bodies{end}.setInitVelocity([0 0 0 0 0 0]');
        end
        
	    model.bodies{end+1} = apbd.BodyRigid(apbd.ShapeCuboid(sides),density);
	    model.bodies{end}.collide = true;
	    model.bodies{end}.mu = mu;
	    %R = se3.aaToMat([1 1 1] / norm([1 1 1]), pi/2);
        R = se3.aaToMat([1 0 0], 0);
	    E = eye(4);
	    x = 0;
	    y = 6*w;
	    z = w;
        E(1:3,1:3) = R;
	    E(1:3,4) = R * [x y z]';
	    model.bodies{end}.setInitTransform(E);
        
        for i = 1:n-1
            if(i<5)
                model.joints{end+1} = JointHinge(model.bodies{i}, model.bodies{i+1}, false, [0 (i-1)*6*w 10*w]' ,[1 0 0]', -10000000*(5-i*2)*ones(model.steps,1));
            else
                model.joints{end+1} = JointHinge(model.bodies{i}, model.bodies{i+1}, false, [0 0 (i+1)*6*w]' ,[1 0 0]');
            end
        end
        groundCollisionList = 1:n;
        groundCollisionList(end+1) = 4;
        bodyCollisionList = [2 4; 3 4;];
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

end
