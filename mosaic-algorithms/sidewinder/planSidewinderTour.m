function [topo_tour, inst_grid, inst_tour, grid_dirx, grid_diry, sweepDir1, sweepDir2] = ...
planSidewinderTour(target, roi, sc, inst, inittime, olapx, olapy, angle)
% This function plans an observation tour using a modified Boustrophedon
% decomposition method. It calculates an optimal path for observing a ROI
% on a target body, considering the spacecraft's starting position.
% We project the ROI's polygon onto the instrument's focal plane. 
% Here, we design the grid based on the image plane's reference. The image
% plane is built according to the FOV's parameters, retrieved from the
% instrument kernel. We also devise the traversal in the image plane (tour), 
% and then transform it back to the topographical coordinates.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [topo_tour, inst_grid, inst_tour, grid_dirx, grid_diry, sweepDir1, sweepDir2] = ...
%                   planSidewinderTour(target, roi, sc, inst, inittime, olapx, olapy, angle)
%
% Inputs:
%   > target:       string name of the target body
%   > roi:          matrix containing the vertices of the uncovered area
%                   of the ROI polygon. The vertex points are expressed in 
%                   2D, in latitudinal coordinates [º]
%   > sc:           string name of the spacecraft
%   > inst:         string name of the instrument
%   > inittime:     start time of the planning horizon, in TDB seconds past
%                   J2000 epoch
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > angle:        observation angle, influencing the orientation of
%                   observation footprints. In other words, it is the angle
%                   that the footprint axes define with respect to the
%                   reference axes of the topographical grid (east-north)
% 
% Outputs:
%   > topo_tour:    tour path in topographical coordinates (lat/lon on the
%                   target body), in [deg]
%   > inst_grid:    grid of potential observation points in instrument
%                   frame coordinates
%   > inst_tour:    tour path in instrument frame coordinates
%   > grid_dirx:    direction of the grid along the x-axis in the
%                   instrument frame
%   > grid_diry:    direction of the grid along the y-axis in the
%                   instrument frame
%   > sweepDir1, sweepDir2: directions defining the sweep of the
%                           Boustrophedon decomposition, derived according 
%                           to the spacecraft's position with respect to
%                           the ROI

% Pre-allocate variables
origin = [0, 0]; % initialize grid origin for grid generation
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2))); % point camera at 
% ROI's centroid

% Project ROI to the instrument plane
targetArea = topo2inst(roi, cx, cy, target, sc, inst, inittime);

% Closest polygon side to the spacecraft's ground track position (this
% will determine the coverage path)
gt1 = groundtrack(sc, inittime, target); % initial ground track position
gt2 = groundtrack(sc, inittime + 500, target); % future ground track position
gt1 = topo2inst(gt1, cx, cy, target, sc, inst, inittime); % projected initial position
gt2 = topo2inst(gt2, cx, cy, target, sc, inst, inittime + 500); % projected future position
% Calculate the closest side of the target area to the spacecraft's ground 
% track, determining the observation sweep direction
[sweepDir1, sweepDir2] = closestSide(gt1, gt2, targetArea, angle);

% Retrieve the field of view (FOV) bounds of the instrument and calculate 
% the dimensions of a reference observation footprint
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
[~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
fpref.width = width;
fpref.height = height;
fpref.angle = angle; % we enforce the orientation angle of the footprint 
% to be that of its projection onto the topographical grid with respect to
% the reference axes (east-north), given as an input. grid2D will use this
% angle to orient the ROI according to this orientation
% [Future work]: This angle could be given as an input to grid2D?

% Focal plane grid discretization based on the reference footprint (FOV
% plane) and specified overlap
[inst_grid, grid_dirx, grid_diry] = grid2D(fpref, olapx, olapy, origin, targetArea);

% Boustrophedon decomposition to generate the grid traversal
inst_tour = boustrophedon(inst_grid, sweepDir1, sweepDir2);

% Convert grid and tour from instrument frame to topographical coordinates
topo_tour = inst2topo(inst_tour, cx, cy, target, sc, inst, inittime);

% Remove empty elements from the tour, which may result from unobservable
% regions within the planned plath
emptyCells = cellfun(@isempty, topo_tour); % find indices of empty cells
topo_tour(emptyCells) = []; % remove empty cells
end