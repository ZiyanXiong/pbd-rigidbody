model = createTestModel(21,1e-2,100);
model.solverType = 1;
model.init();
model.simulate();


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

