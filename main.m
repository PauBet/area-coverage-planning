clc; clear all; close all;
% Paula Betriu - August 2022
% Approximation Planning Algorithms
% This program implements a series of area coverage algorithms, based on 
% [1]. These algorithms aim to build the optimal plan for an instrument
% to represent or cover a specific area from a body surface.
%
% [1] Shao E., Byon A., Davies C., Davis E., Knight R., Lewellen G., 
% Trowbridge M. and Chien S. Area Coverage Planning with 3-axis Steerable,
% 2D Framing Sensors. 2018.

% Load mission info (kernels, SPICE ids, etc.)
input_data;

% Pre-allocation of variables... 
%inittime = cspice_str2et('1997 DEC 16 10:00:00.000 TDB'); % mosaic start
stoptime = cspice_str2et('1998 MAY 30 00:00:00.000 TDB'); % mosaic end (max)

% Regions of interest
% Annwn Regio [lon, lat] = [20, 40]º
roi{1} = [50 35;
          50  5;
          30  5;
          30 35;];

inittime{1} = cspice_str2et('1998 MAR 29 12:44:12.000 TDB'); % closest approach
%inittime = cspice_str2et('1997 DEC 16 10:30:00.000 TDB'); % mosaic start
    
% Pwyll crater [lon, lat] = [90, -25]º;
%ax = mapPlot('europa-map-shifted.png');
%roi{2} = drawpolygon(ax);
roi{2} = [ 78   -18;
           90   -14;
          109   -13;
          101   -23;
           99   -30;
           85   -33;
           71   -28;];

inittime{2} = cspice_str2et('1998 MAR 29 12:38:02.000 TDB'); % closest approach

% Tara Regio [lon, lat] = [-75, -10]º;
roi{3} = [-55   0;
          -85   0;
          -85 -20;
          -55 -20;]; % roi of roi polygon

inittime{3} = cspice_str2et('1998 MAR 29 14:20:52.000 TDB'); % closest approach

% Cilix crater [lon, lat] = [180, 0]º;
roi{4} = [-177  3;
          -177 -3;
           177 -3;
           177  3;];

inittime{4} = cspice_str2et('1998 MAR 29 13:30:22.000 TDB'); % closest approach

tobs  = 10; % [s] between observations
olapx = 10; % [%] of overlap in x direction
olapy = 10; % [%] of overlap in y direction
videosave = 0; % =1 saves a video with the evolution of the coverage map
% along the observation plan

% Coverage figure:
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex planetary bodies
ax = mapPlot('europa-map-shifted.png');
c  = {'g', 'y', 'c', 'm'};
for i=1:length(roi)
    [x, y] = amsplit(roi{i}(:, 1), roi{i}(:, 2));
    plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
        c{i}, 'linewidth', 1)
end
legend('Annwn Regio', 'Pwyll Crater', 'Tara Regio', 'Footprint',...
    'Ground track', 'AutoUpdate', 'off')

%% Mosaic algorithms

for i=2:2
    
    xlim([min(roi{i}(:, 1))-10  max(roi{i}(:, 1)) + 10])
    ylim([min(roi{i}(:, 2))-10  max(roi{i}(:, 2)) + 10])

    % Sidewinder
%    [A, cv] = sidewinder(inittime{i}, stoptime, tobs, inst, sc, target, ...
%       roi{i}, olapx, olapy, ax, c{i}, videosave);

    % Replanning Sidewinder
     [A, cv] = replanningSidewinder(inittime{i}, stoptime, tobs, inst, ...
         sc, target, roi{i}, olapx, olapy, ax, c{i}, 0);

    % Online Frontier Repair

    % Grid Nibbler
    % A = gridNibbler(startTime, endTime, step, instName, scName, targetName, roi, 20, 20, 1);

    % Re-plot the ROI (for aesthetic purposes)
    [x, y] = amsplit(roi{i}(:, 1), roi{i}(:, 2));
    plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
        c{i}, 'linewidth', 1)

    % Print coverage
    fprintf('Total coverage: %.3f%% \n', cv*100);

end
xlim([-180 180])
ylim([-90 90])

% End SPICE
endSPICE;