clc; close all; clear all;
% Revision of Sidewinder:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see
% Re-implementation of sidewinder function with these new feature (and
% other code improvements)

% Load mission info (kernels, SPICE ids, etc.)
input_data_fig1;

% Sidewinder
[A, fpList] = sidewinder(roistruct(1).inittime, ...
    stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3*1e-3, 0);

% Plot tour
plotTour(A, fpList, roistruct, sc, target)

% FOM post-process
post_process_fig1;

% Save figure [PDF]
figpath = '.';
set(gcf, 'Units', 'inches', 'Position', [3,3,9.2,6]);
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
filename = fullfile(figpath, 'fig1_sidewinder');
print(gcf, filename, '-dpdf', '-r1200')

% % Zoom in
% xlim([120 180])
% ylim([0 40])
% xtick = 120:15:180;
% ytick = 0:10:40;
% 
% % x tick label
% xtickstr = string();
% for i=1:length(xtick)
%     if xtick(i) < 0 && xtick(i) > -180
%         xtickstr(i) = strcat(num2str(-xtick(i)), 'º', 'W');
%     elseif xtick(i) > 0 && xtick(i) < 180
%         xtickstr(i) = strcat(num2str(xtick(i)), 'º', 'E');
%     else
%         xtickstr(i) = strcat(num2str(abs(xtick(i))), 'º');
%     end
% end
% 
% % y tick label
% ytickstr = string();
% for i=1:length(ytick)
%     if ytick(i) < 0
%         ytickstr(i) = strcat(num2str(-ytick(i)), 'º', 'S');
%     elseif ytick(i) > 0
%         ytickstr(i) = strcat(num2str(ytick(i)), 'º', 'N');
%     else
%         ytickstr(i) = strcat(num2str(ytick(i)), 'º');
%     end
% end
% set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
%     'YTickLabel', ytickstr, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'Layer', 'top', 'GridColor', 'w', 'GridAlpha', 1, 'GridLineStyle', ':');
% 
% % Save figure [PDF]
% figpath = '.';
% set(gcf, 'Units', 'inches', 'Position', [3,3,7.4,5.3]);
% %set(gcf, 'Units', 'inches', 'Position', [3,3,9.2,6]);
% pos = get(gcf, 'Position');
% set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
%     'PaperSize', [pos(3), pos(4)]);
% filename = fullfile(figpath, 'sidewinder');
% print(gcf, filename, '-dpdf', '-r600')