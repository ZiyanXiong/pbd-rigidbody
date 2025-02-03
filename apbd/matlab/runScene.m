rng(42);
%{
model = createTestModel(4,1e-4,1);
model.init();
model.simulate();
%}

for i = 13
    model = SIG25Model(i,1/100,150,2);
    model.init();
    model.drawHz = 0;
    model.simulate();
    if(model.solverType == 1)
        fileName = "iterVec_rVec_TGS.mat";
    elseif(model.solverType == 2)
        fileName = "iterVec_rVec_2PSP.mat";
    else
        fileName = sprintf("iterVec_rVec_TGS_%d.mat", model.substeps);
    end
    
    iterVec = model.iterVec;
    rVec = model.rVec;
    save(strcat(model.resultFolder, fileName), 'iterVec', 'rVec');
    %SIG25Plot(i, false);
end


%{
for i = 1
    model = SIG25Model(i,1/60,150,2);
    model.drawHz = 60;
    model.init();
    model.simulate();
    if(model.solverType == 1)
        fileName = sprintf("iterVec_rVec_TGS_%d.mat", model.substeps);
    elseif(model.solverType == 2)
        fileName = "iterVec_rVec_2PSP.mat";
    else
        fileName = sprintf("iterVec_rVec_TGS_%d.mat", model.substeps);
    end
    
    iterVec = model.iterVec;
    rVec = model.rVec;
    save(strcat(model.resultFolder, fileName), 'iterVec', 'rVec');

    for j = [50,150,500]
        model = SIG25Model(i,1/60,j,3);
        model.drawHz = 1;
        model.init();
        model.simulate();
        if(model.solverType == 1)
            fileName = sprintf("iterVec_rVec_TGS_%d.mat", model.substeps);
        elseif(model.solverType == 2)
            fileName = "iterVec_rVec_2PSP.mat";
        else
            fileName = sprintf("iterVec_rVec_TGS_%d.mat", model.substeps);
        end
        
        iterVec = model.iterVec;
        rVec = model.rVec;
        save(strcat(model.resultFolder, fileName), 'iterVec', 'rVec');
    end
    fprintf("Finished scene %d \n", i);
end
%}

%{
resultFolder = strcat('Results\MATLAB\');
if ~exist(resultFolder, 'dir')
   mkdir(resultFolder)
end
model.video = VideoWriter(strcat(resultFolder, ...
                'Solver_GS', ...
                '_timestep_', string(model.h), '.avi'));
open(model.video);
model.init();
model.simulate();
close(model.video);
%}

%% Code for generating video results
%{
for i = 10:10
    filename = strcat('Scene_',string(i));
    resultFolder = strcat('Results\MATLAB\');
    if ~exist(resultFolder, 'dir')
       mkdir(resultFolder)
    end
    for solvertype = 3:-1:2
        for substeps = [20 50 100]
            model = createTestModel(i, 1e-2, substeps);
            model.solverType = solvertype;
            model.savedBodyStatesPath = strcat(resultFolder, filename, ...
                           '_timestep_', string(model.h), ...
                            '_substeps_', string(model.substeps), '.txt');
            model.init();
            model.simulate();
        end
    end
end
%}

%{
for i = 13:14
    filename = strcat('Scene',string(i),'\');
    resultFolder = strcat('Results\',filename);
    if ~exist(resultFolder, 'dir')
       mkdir(resultFolder)
    end
    for solvertype = 1:3
        for substeps = [20 50 100]
            model = createTestModel(i, 1e-2, substeps);
            model.solverType = solvertype;
            model.video = VideoWriter(strcat(resultFolder, ...
                            'Solver_', string(model.solverType), ...
                            '_timestep_', string(model.h), ...
                            '_substeps_', string(model.substeps), '.avi'));
            open(model.video);
            model.init();
            model.simulate();
            close(model.video);
        end
    end
end
%}

