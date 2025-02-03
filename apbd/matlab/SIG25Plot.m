function SIG25Plot(modelID, animation)
    resultFolder = sprintf("Results\\Scene\\%d\\", modelID);
    plotFolder = strcat(resultFolder, "plots\\");
    if ~exist(plotFolder, 'dir')
       mkdir(plotFolder)
    end

    %{
    fileName = "iterVec_rVec_GS.mat";
    load(strcat(resultFolder, fileName), "rVec", "iterVec");
    GS_rVec = rVec';
    GS_iterVec = iterVec';

    fileName = "iterVec_rVec_2PSP.mat";
    load(strcat(resultFolder, fileName), "rVec", "iterVec");
    TPSP_rVec = rVec';
    TPSP_iterVec = iterVec';
    
    for i = 1:steps
        load(strcat(resultFolder,sprintf("\\residual_per_iteration\\2PSP_%d.mat",i-1)), "rs_gs", "rs_2psp");
        f = figure('visible','off');
        semilogy(1:size(rs_gs,1), rs_gs, '-', 'DisplayName','Gauss-Seidel','linewidth',1.5);
        hold on;
        semilogy(1:size(rs_2psp,1), rs_2psp,'-', 'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
        ylabel('Residual');    % Label for left y-axis
        xlabel('Iteration number');           % Label for x-axis
        legend;
        title(sprintf('Residul after each iteration at time step %d', i));
        grid on;
        exportgraphics(f, strcat(plotFolder, sprintf("residual_per_iteration_%d.tif", i)), 'Resolution',110);
        close(f);

        f = figure('visible','off');
        semilogy(1:size(GS_rVec,1), GS_rVec, '-', 'DisplayName','Gauss-Seidel','linewidth',1.5);
        hold on;
        semilogy(1:size(TPSP_rVec,1), TPSP_rVec, '-' ,'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
        % Add vertical red line at x = 10
        xline(i, '--k','HandleVisibility', 'off', 'LineWidth', 1.5);
        ylabel('Residual after solve');    % Label for left y-axis
        xlabel('Time step');           % Label for x-axis
        title('Residul at the end of each time step');
        legend();
        grid on;
        exportgraphics(f, strcat(plotFolder, sprintf("residual_per_timestep_%d.tif", i)), 'Resolution',110);
        close(f);

        f = figure('visible','off');
        plot(1:size(GS_iterVec,1), GS_iterVec, '-', 'DisplayName','Gauss-Seidel','linewidth',1.5);
        hold on;
        plot(1:size(TPSP_iterVec,1), TPSP_iterVec, '-', 'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
        xline(i, '--k','HandleVisibility', 'off', 'LineWidth', 1.5);
        ylim([0 1500])
        ylabel('Total Iteration number');    % Label for left y-axis
        xlabel('Time step');           % Label for x-axis
        title('Iteration number of each time step');
        legend;
        grid on;
        exportgraphics(f, strcat(plotFolder, sprintf("iteration_per_timestep_%d.tif", i)), 'Resolution',110);
        close(f);
    end
    %}

    load(strcat(resultFolder,sprintf("\\residual_per_iteration\\2PSP_1.mat")), "rs_gs", "rs_2psp");
    f = figure('visible','off');
    semilogy(1:size(rs_gs,1), rs_gs, '-', 'DisplayName','Gauss-Seidel','linewidth',1.5);
    hold on;
    semilogy(1:size(rs_2psp,1), rs_2psp,'-', 'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
    ylabel('Residual');    % Label for left y-axis
    xlabel('Iteration number');           % Label for x-axis
    legend;
    title(sprintf('Residul after each iteration at time step 1'));
    grid on;
    exportgraphics(f, strcat(plotFolder, sprintf("residual_per_iteration_1.tif")), 'Resolution',200);
    savefig(f, strcat(plotFolder, sprintf("residual_per_iteration_1.fig")));
    close(f);


    fileName = "iterVec_rVec_TGS_50.mat";
    load(strcat(resultFolder, fileName), "rVec", "iterVec");
    TGS_50_rVec = rVec';
    TGS_50_iterVec = iterVec';

    fileName = "iterVec_rVec_TGS_150.mat";
    load(strcat(resultFolder, fileName), "rVec", "iterVec");
    TGS_150_rVec = rVec';
    TGS_150_iterVec = iterVec';

    fileName = "iterVec_rVec_TGS_500.mat";
    load(strcat(resultFolder, fileName), "rVec", "iterVec");
    TGS_500_rVec = rVec';
    TGS_500_iterVec = iterVec';

    fileName = "iterVec_rVec_2PSP.mat";
    load(strcat(resultFolder, fileName), "rVec", "iterVec");
    TPSP_rVec = rVec';
    TPSP_iterVec = iterVec';

    steps = length(rVec);


    if(animation)
        for i = 1:steps
            f = figure('visible','off');
            semilogy(1:size(TGS_50_rVec,1), TGS_50_rVec, '-', 'DisplayName','Gauss-Seidel with 50 sub steps','linewidth',1.5);
            hold on;
            semilogy(1:size(TGS_150_rVec,1), TGS_150_rVec, '-' ,'DisplayName','Gauss-Seidel with 150 sub steps','linewidth',1.5);
            hold on;
            semilogy(1:size(TGS_500_rVec,1), TGS_500_rVec, '-' ,'DisplayName','Gauss-Seidel with 500 sub steps','linewidth',1.5);
            hold on;
            semilogy(1:size(TPSP_rVec,1), TPSP_rVec, '-' ,'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
    
            % Add vertical red line at x = 10
            xline(i, '--k','HandleVisibility', 'off', 'LineWidth', 1.5);
            ylabel('Residual after solve');    % Label for left y-axis
            xlabel('Time step');           % Label for x-axis
            title('Residul at the end of each time step');
            legend('Location', 'best');
            grid on;
            exportgraphics(f, strcat(plotFolder, sprintf("residual_per_timestep_%d.tif", i)), 'Resolution',200);
            close(f);
    
            f = figure('visible','off');
            plot(1:size(TGS_50_iterVec,1), TGS_50_iterVec, '-', 'DisplayName','Gauss-Seidel with 50 sub steps','linewidth',1.5);
            hold on;
            plot(1:size(TGS_150_iterVec,1), TGS_150_iterVec, '-', 'DisplayName','Gauss-Seidel with 150 sub steps','linewidth',1.5);
            hold on;
            plot(1:size(TGS_500_iterVec,1), TGS_500_iterVec, '-', 'DisplayName','Gauss-Seidel with 500 sub steps','linewidth',1.5);
            hold on;
            plot(1:size(TPSP_iterVec,1), TPSP_iterVec, '-', 'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
            xline(i, '--k','HandleVisibility', 'off', 'LineWidth', 1.5);
            ylim([0 1000])
            ylabel('Total Iteration number');    % Label for left y-axis
            xlabel('Time step');           % Label for x-axis
            title('Iteration number / Sub steps of each time step');
            legend;
            grid on;
            exportgraphics(f, strcat(plotFolder, sprintf("iteration_per_timestep_%d.tif", i)), 'Resolution',200);
            close(f);
        end
    end

        f = figure('visible','off');
        semilogy(1:size(TGS_50_rVec,1), TGS_50_rVec, '-', 'DisplayName','Gauss-Seidel with 50 sub steps','linewidth',1.5);
        hold on;
        semilogy(1:size(TGS_150_rVec,1), TGS_150_rVec, '-' ,'DisplayName','Gauss-Seidel with 150 sub steps','linewidth',1.5);
        hold on;
        semilogy(1:size(TGS_500_rVec,1), TGS_500_rVec, '-' ,'DisplayName','Gauss-Seidel with 500 sub steps','linewidth',1.5);
        hold on;
        semilogy(1:size(TPSP_rVec,1), TPSP_rVec, '-' ,'DisplayName','2 Pass Shock Propagation','linewidth',1.5);

        % Add vertical red line at x = 10
        ylabel('Residual after solve');    % Label for left y-axis
        xlabel('Time step');           % Label for x-axis
        title('Residul at the end of each time step');
        legend('Location', 'best');
        grid on;
        exportgraphics(f, strcat(plotFolder, sprintf("residual_per_timestep.tif")), 'Resolution',200);
        savefig(f, strcat(plotFolder, sprintf("residual_per_timestep.fig")));
        close(f);

        f = figure('visible','off');
        plot(1:size(TGS_50_iterVec,1), TGS_50_iterVec, '-', 'DisplayName','Gauss-Seidel with 50 sub steps','linewidth',1.5);
        hold on;
        plot(1:size(TGS_150_iterVec,1), TGS_150_iterVec, '-', 'DisplayName','Gauss-Seidel with 150 sub steps','linewidth',1.5);
        hold on;
        plot(1:size(TGS_500_iterVec,1), TGS_500_iterVec, '-', 'DisplayName','Gauss-Seidel with 500 sub steps','linewidth',1.5);
        hold on;
        plot(1:size(TPSP_iterVec,1), TPSP_iterVec, '-', 'DisplayName','2 Pass Shock Propagation','linewidth',1.5);
        ylim([0 1000])
        ylabel('Total Iteration number');    % Label for left y-axis
        xlabel('Time step');           % Label for x-axis
        title('Iteration number / Sub steps of each time step');
        legend;
        grid on;
        exportgraphics(f, strcat(plotFolder, sprintf("iteration_per_timestep.tif")), 'Resolution',200);
        savefig(f, strcat(plotFolder, sprintf("iteration_per_timestep.fig")));
        close(f);
end