% Representative examples of the approximation heuristics performance
% July 2023
clc; clear all; close all;

% Load mission info (kernels, SPICE ids, etc.)
inputdata;

% Coverage figure:
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex planetary bodies
ax = mapPlot('ganymede-map.jpg');

%% Mosaic algorithms
v2 = [];
obsDic = dictionary();
%     if abs(min(roi{i}(:, 1)) -  max(roi{i}(:, 1))) > 180
%         xlim([175  180])
%         ylim([min(roi{i}(:, 2)) - 10  max(roi{i}(:, 2)) + 10])
%     else
%         xlim([min(roi{i}(:, 1)) - 10  max(roi{i}(:, 1)) + 10])
%         ylim([min(roi{i}(:, 2)) - 10  max(roi{i}(:, 2)) + 10])
%     end

% Pre-allocate variables
slew_rate = 3*1e-3; % [rad/s]
t = zeros(length(roistruct), 1);
coverage = zeros(length(roistruct), 1);
overlap = zeros(length(roistruct), 1);
makespan = zeros(length(roistruct), 1);
strname = string();
nfp = zeros(length(roistruct), 1);

% Mosaic algorithms
for i=1:length(roistruct)
    % Get roi polygon and initial time observation...
    roi = roistruct(i).vertices;
    inittime = roistruct(i).inittime;
    strname(i) = roistruct(i).name;

    [A, fplist] = frontierRepair(inittime, stoptime, ...
        tobs, inst, sc, target, roi, olapx, olapy, slew_rate, 'highres');

    % Get coverage, overlap and makespan
    [coverage(i), overlap(i)] = roicoverage(target, roi, fplist);
    makespan(i) = fplist(end).t + tobs - fplist(1).t;
    nfp(i) = length(fplist);

    % Plot tour
    plotTour_m(A, fplist, roistruct, target, ax);
    drawnow
end
title(tilealg + " coverage map")
% Save figure [PDF]
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
filename = [tilealg, '_coverage_map'];
print(gcf, filename, '-dpdf', '-r600')

varnames = ["ROI", "Coverage", "Overlap", "Num Footprints", "Makespan", "CPU Time"];
tab = table(strname', coverage, overlap, nfp, makespan, t, 'VariableNames', varnames);
disp(tab)

% End SPICE
endSPICE;