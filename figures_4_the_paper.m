% Representative examples of the approximation heuristics performance
% July 2023
clc; clear all; close all

% Load mission info (kernels, SPICE ids, etc.)
input_data;

% Choose mosaic algorithm: 'sidewinder', 'r_sidewinder', 'onlinefrontier',
% 'gridnibbler'
tilealg = 'sidewinder';

% Pre-allocation of variables... 
stoptime = cspice_str2et('1998 MAY 30 00:00:00.000 TDB'); % mosaic end (max)
tobs  = 10; % [s] between observations
olapx = 20; % [%] of overlap in x direction
olapy = 20; % [%] of overlap in y direction
videosave = 0; % =1 saves a video with the evolution of the coverage map
% along the observation plan

% Regions of interest
% Annwn Regio [lon, lat] = [40, 20]º
roi = [55 30;
       55 10;
       35 10;
       35 30;];

inittime = cspice_str2et('1998 MAR 29 12:53:00.000 TDB'); % closest approach
region = 'Annwn Regio';

% Coverage figure:
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex planetary bodies
ax = mapPlot('europa-map.jpg');
c  = {'g', 'y', 'c', 'm', 'r', 'k'};
% Define the RdYlBu colormap
map = getPyPlot_cMap('tab20', 7, [], '/usr/local/bin/python3');
[x, y] = amsplit(roi(:, 1), roi(:, 2));
plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
    [0.09,0.09,0.54], 'linewidth', 2)

%% Mosaic algorithms
if videosave
    v2 = VideoWriter('topography_map_sidewinder', 'MPEG-4');
    v2.FrameRate = 2;
    open(v2)
    writeVideo(v2, getframe(gcf));
else
    v2 = [];
end
obsDic = dictionary();
%     if abs(min(roi{i}(:, 1)) -  max(roi{i}(:, 1))) > 180
%         xlim([175  180])
%         ylim([min(roi{i}(:, 2)) - 10  max(roi{i}(:, 2)) + 10])
%     else
%         xlim([min(roi{i}(:, 1)) - 10  max(roi{i}(:, 1)) + 10])
%         ylim([min(roi{i}(:, 2)) - 10  max(roi{i}(:, 2)) + 10])
%     end

switch tilealg
    case 'sidewinder'
        [A, cv, fplist] = sidewinder(inittime, stoptime, tobs, ...
            inst, sc, target, roi, olapx, olapy, ax, v2, 0);

    case 'r_sidewinder'
        [A, cv, fplist] = replanningSidewinder(inittime, ...
            stoptime, tobs, inst, sc, target, roi, olapx, olapy, ...
            ax, []);

    case 'onlinefrontier'
        [A, cv, fplist] = frontierRepair(inittime, stoptime, ...
            tobs, inst, sc, target, roi, olapx, olapy, ax, ...
            'g', v2);

    case 'gridnibbler'
        A = gridNibbler(startTime, endTime, step, instName, scName, ...
            targetName, roi, 20, 20, 1);

    case 'neighbours'
        A = neighbour_placement(startTime, tobs, inst, '-77', ...
            target, roi, ax);
end

clear grid2D planSidewinderTour;

% Re-plot the ROI (for aesthetic purposes)
[x, y] = amsplit(roi(:, 1), roi(:, 2));
plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
    [0.09,0.09,0.54], 'linewidth', 1)
set(gca, 'xlim', [-180 180], 'ylim', [-90 90])
legend('Region of interest', 'Footprint', 'Start point', 'Ground track', 'AutoUpdate', 'off', 'location', 'northwest')
if videosave
    writeVideo(v2, getframe(gcf));
    close(v2);
end

% Zoom in
xlim([15 75])
ylim([-5 35])
xtick = 15:15:75;
ytick = -5:10:35;

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
    'YTickLabel', ytickstr, 'Layer', 'top', 'GridColor', 'w', ...
    'GridAlpha', 1, 'GridLineStyle', ':');

% End SPICE
endSPICE;