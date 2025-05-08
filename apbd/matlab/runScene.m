for i = 3
    model = JointModel(i,1/100,10,1);
    model.init();
    model.drawHz = 100;
    model.simulate();
    if(model.solverType == 1)
        fileName = sprintf("iterVec_rVec_TGS_%d.mat", model.substeps);
    elseif(model.solverType == 2)
        fileName = "iterVec_rVec_GPQP.mat";
    end
    
    iterVec = model.iterVec;
    rVec = model.rVec;
    save(strcat(model.resultFolder, fileName), 'iterVec', 'rVec'); 
    %SIG25Plot(i, false);
end

