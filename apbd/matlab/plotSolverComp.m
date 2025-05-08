for k = 1:1
     processPlot(k,10);
end

%%
function processPlot(scene, step)
fprintf('Process Scene %d\n',scene);
resultFolder = sprintf("Results\\Scene\\%d\\", scene);
fileName = sprintf("rVec_GPQP_step_%d.mat", step);
load(strcat(resultFolder, fileName), "rs_gpqp");

fileName = sprintf("rVec_GS_step_%d.mat", step);
load(strcat(resultFolder, fileName), "rs_gs");


%f = figure('visible','off');
f = figure();
semilogy(1:size(rs_gpqp,1), rs_gpqp, '-', 'DisplayName','GPQP','linewidth',2);
hold on;
semilogy(1:size(rs_gs,1), rs_gs, '-', 'DisplayName','GS','linewidth',2);
hold on;
legend('Location', 'bestoutside');
set(gca,'FontSize',22);
title('');
ylabel('Residual');
xlabel('Iterations');
ylim([1e-4,1e3]);
ax = gca; % Get the current axes
ax.XTick = 0:100:size(rs_gs,1); % Major ticks
ax.YTick = logspace(-6,4,6); % Major ticks
grid on;
exportgraphics(f, strcat(resultFolder,sprintf("residuals_plot_scene_%d_step_%d.png",scene,step)), 'Resolution',400);
savefig(f, strcat(resultFolder,sprintf("residuals_plot_scene_%d_step_%d.fig",scene,step)));
close(f);

end

