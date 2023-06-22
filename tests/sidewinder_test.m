clc; clear all; close all
% Paula Betriu - April 2023
% Tiling Algorithms test of robustness

% Load mission info (kernels, SPICE ids, etc.)
input_data;

% Choose mosaic algorithm: 'sidewinder', 'r_sidewinder', 'onlinefrontier',
% 'gridnibbler'
tilealg = 'sidewinder';

% Get the region of interest
figure
hold on; grid minor; axis equal;
set(gca, 'xlim', [-180 180], 'ylim', [-90 90], 'fontsize', 20)
ax = gca;
roipoly = drawpolygon(ax);
roi = roipoly.Position;

% Pre-allocation of variables... 
tobs  = 10; % [s] between observations
olapx = 10; % [%] of overlap in x direction
olapy = 10; % [%] of overlap in y direction

% Test time window
wstart = cspice_str2et('1998 MAR 29 11:30:00.000 TDB');
wend   = cspice_str2et('1998 MAR 29 12:30:00.000 TDB');
step   = 3*60;
et     = wstart:step:wend;

% Coverage figure:
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex planetary bodies
for i=1:length(et)
    ax = mapPlot('europa-map-shifted.png');
    plot(ax, polyshape(roi(:, 1), roi(:, 2)), 'FaceColor', 'none', 'EdgeColor', 'w', ...
        'linewidth', 2)
    xlim([min(roi(:, 1)) - 10  max(roi(:, 1)) + 10])
    ylim([min(roi(:, 2)) - 10  max(roi(:, 2)) + 10])
    [A, cv, fplist] = sidewinder(et(i), wend + 4000, tobs, ...
        inst, sc, target, roi, olapx, olapy, ax, 'b', []);
    clear grid2D planSidewinderTour;
    delete(gcf)
end

% End SPICE
endSPICE;