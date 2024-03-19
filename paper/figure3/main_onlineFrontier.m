clc; close all; clear all;
% Revision of grid discretization:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see
% Re-implementation of sidewinder function with these new feature (and
% other code improvements)

% Load mission info (kernels, SPICE ids, etc.)
input_data_fig3;
mosaic = 'onlinefrontier';
roiname = char(lower(roistruct(1).name));
roiname(isspace(roiname)) = [];
name = ['post_process_',roiname];

% Online Frontier
[A, fpList] = frontierRepair(roistruct(1).inittime, ...
    stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3*1e-3);

% Plot tour
plotTour(A, fpList, roistruct, sc, target)
title(roistruct(1).name)
leg = legend('NumColumns', 2, 'Location', 'north');
leg.String(end) = [];

% FOM post-process
run(name);