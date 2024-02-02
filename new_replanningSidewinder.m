clc; close all; clear all;
% Revision of Sidewinder:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see
% Re-implementation of sidewinder function with these new feature (and
% other code improvements)

% % Load mission info (kernels, SPICE ids, etc.)
% input_data;
% %inittime = cspice_str2et('1998 MAR 29 12:53:00.000 TDB'); % closest approach
% inittime = cspice_str2et('1998 MAR 29 12:55:00.000 TDB'); % closest approach
% stoptime = inittime + 24*3600;
% tcadence = 10; % observation cadence [s]
% olapx = 20; olapy = 20; % overlap in [%]
% speedUp = 0;
% 
% % Define ROI (Annwn Regio)
% roi = [65 30;
%        65 15;
%        45 15;
%        45 30;];

% Load mission info (kernels, SPICE ids, etc.)
input_data;
%inittime = cspice_str2et('1998 MAR 29 12:53:00.000 TDB'); % closest approach
inittime = cspice_str2et('1998 MAR 29 12:35:00.000 TDB'); % closest approach
stoptime = inittime + 24*3600;
tcadence = 20; % observation cadence [s]
olapx = 20; olapy = 20; % overlap in [%]
speedUp = 0;

% Define ROI (Annwn Regio)
roi = [25  25;
       25  -5;
       -5  -5;
       -5  25;];

% Point camera at ROI's centroid
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(1).vertices = roi;
roistruct(1).cpoint = [cx, cy];
roistruct(1).name = "Annwn Regio";

% Sidewinder
tic
[A, fpList] = replanningSidewinder2(inittime, ...
    stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3*1e-3, 0);
toc
% Plot tour
plotTour(A, fpList, roistruct, sc, target)

% Zoom in
xlim([-30 60])
ylim([-15 35])
xtick = -30:15:60;
ytick = -15:10:35;

% x tick label
xtickstr = string();
for i=1:length(xtick)
    if xtick(i) < 0 && xtick(i) > -180
        xtickstr(i) = strcat(num2str(-xtick(i)), 'º', 'W');
    elseif xtick(i) > 0 && xtick(i) < 180
        xtickstr(i) = strcat(num2str(xtick(i)), 'º', 'E');
    else
        xtickstr(i) = strcat(num2str(abs(xtick(i))), 'º');
    end
end

% y tick label
ytickstr = string();
for i=1:length(ytick)
    if ytick(i) < 0
        ytickstr(i) = strcat(num2str(-ytick(i)), 'º', 'S');
    elseif ytick(i) > 0
        ytickstr(i) = strcat(num2str(ytick(i)), 'º', 'N');
    else
        ytickstr(i) = strcat(num2str(ytick(i)), 'º');
    end
end
set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
    'YTickLabel', ytickstr, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Layer', 'top', 'GridColor', 'w', 'GridAlpha', 1, 'GridLineStyle', ':');

% Save figure [PDF]
figpath = '.';
set(gcf, 'Units', 'inches', 'Position', [3,3,9.2,6]);
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
filename = fullfile(figpath, 'r_sidewinder');
print(gcf, filename, '-dpdf', '-r1200')