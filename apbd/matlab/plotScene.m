for k = 2:12
     processIter(k);
end

for k = 2:12
    processRes(k);
end

%%
function processIter(scene)
fprintf('Process Scene %d\n',scene);
resultFolder = sprintf("Results\\Scene\\%d\\", scene);
plotFolder = strcat(resultFolder, "plots\\");
for i = 1: 60
    openfig(strcat(plotFolder,'iteration_per_timestep.fig'),'invisible');
    lines = findall(gca, 'Type', 'Line');
    set(lines, 'LineWidth', 2);
    legend('GS-50','GS-150','GS-500','2PSP', 'Location', 'bestoutside');
    set(gca,'FontSize',22);
    xline(i, '--k','HandleVisibility', 'off', 'LineWidth', 1.5);
    title('');
    ylabel('Iters');
    xlabel('Time Step');
    ylim([0,550]);
    ax = gca; % Get the current axes
    ax.XTick = 0:10:60; % Major ticks
    ax.YTick = 0:100:600; % Major ticks
    set(gcf, 'PaperUnits', 'inches'); % Set units to inches
    set(gcf, 'PaperPosition', [0, 0, 5, 4]); % [left, bottom, width, height]
    %saveas(gcf,['Scene_',scene,'_Iter.pdf']);
    exportgraphics(gcf, strcat(plotFolder, sprintf("iteration_per_timestep_%d.tif", i)), 'Resolution',300);
end
end

%%
function processRes(scene)
fprintf('Process Scene %d\n',scene);
resultFolder = sprintf("Results\\Scene\\%d\\", scene);
plotFolder = strcat(resultFolder, "plots\\");
for i = 1:60
    openfig(strcat(plotFolder,'residual_per_timestep.fig'),'invisible');
    lines = findall(gca, 'Type', 'Line');
    set(lines, 'LineWidth', 2);
    xline(i, '--k','HandleVisibility', 'off', 'LineWidth', 1.5);
    legend off;
    set(gca,'FontSize',22);
    title('');
    ylabel('Residual');
    xlabel('Time Step');
    ylim([1e-5,1e4]);
    ax = gca; % Get the current axes
    ax.XTick = 0:10:60; % Major ticks
    ax.YTick = logspace(-6,4,6); % Major ticks
    set(gcf, 'PaperUnits', 'inches'); % Set units to inches
    set(gcf, 'PaperPosition', [0, 0, 5, 4]); % [left, bottom, width, height]
    %saveas(gcf,['Scene_',scene,'_Res.pdf']);
    exportgraphics(gcf, strcat(plotFolder, sprintf("residual_per_timestep_%d.tif", i)), 'Resolution',300);
end
end
