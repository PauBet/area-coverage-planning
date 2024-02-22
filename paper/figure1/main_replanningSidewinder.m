clc; close all; clear all;
% Revision of grid discretization:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see
% Re-implementation of sidewinder function with these new feature (and
% other code improvements)

% Load mission info (kernels, SPICE ids, etc.)
input_data_fig1;

% Sidewinder
[A, fpList] = replanningSidewinder2(roistruct(1).inittime, ...
    stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3*1e-3, speedUp);
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
filename = fullfile(figpath, 'fig1_replansidewinder');
print(gcf, filename, '-dpdf', '-r1200')