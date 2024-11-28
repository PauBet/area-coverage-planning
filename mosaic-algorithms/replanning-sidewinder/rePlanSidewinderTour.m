function [topo_tour, topo_grid, emptySeed] = rePlanSidewinderTour(target, roi, sc, inst, et, ...
    olapx, olapy, angle, cx, cy, topo_tour, topo_grid, origin_topo, old_origin_topo)
% This function updates the planning of an observation tour by 
% recalculating the grid and the sequence of observation points (tour) 
% based on the current state of the region of interest, instrument 
% orientation, and other parameters. It ensures that the tour remains 
% efficient and covers the ROI effectively, even as conditions change.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        
%
% Inputs:
%   > target:       string name of the target body
%   > roi:          matrix containing the vertices of the uncovered area
%                   of the ROI polygon. The vertex points are expressed in 
%                   2D, in latitudinal coordinates [º]
%   > sc:           string name of the spacecraft
%   > inst:         string name of the instrument
%   > et:           observation time, in TDB seconds past J2000 epoch
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > angle:        observation angle, influencing the orientation of
%                   observation footprints. In other words, it is the angle
%                   that the footprint axes define with respect to the
%                   reference axes of the topographical grid (east-north)
%   > cx, cy:       centroid coordinates of the original roi, used for 
%                   alignment and reference
%   > topo_grid:    cell array of the stare points in topographic 
%                   coordinates. Each cell contains a 2D point [lon, lat], 
%                   in [deg]
%   > origin_topo:  current origin of the topographical grid in latitudinal 
%                   coordinates
%   > old_origin_topo: previous origin of the topographical grid
% 
% Outputs:
%   > topo_tour:    tour path in topographical coordinates (lat/lon on the
%                   target body), in [deg]
%   > topo_grid:    updated grid after replanning
%   > emptySeed:    boolean flag indicating if the seed point (origin_topo)
%                   resulted in an unproductive update, leading to a 
%                   sterile or invalid observation point requiring removal
%                   from the tour

% Pre-allocate variables
persistent fpref; % reference footprint for observations. It corresponds to
% the FOV plane
persistent sweepDir1;
persistent sweepDir2;
if isempty(origin_topo)
    % Point camera at ROI's centroid
    [origin_topo(1), origin_topo(2)] = centroid(polyshape(roi(:, 1), roi(:, 2)));
end

% Build reference tile (it's always going to be the same in the subsequent
% calls)
if isempty(fpref)
    [~, ~, ~, bounds] = ...
        cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
    [~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
    fpref.width = width;
    fpref.height = height;
    fpref.angle = angle;
end

% Intersect ROI with focal plane
targetArea = topo2inst(roi, cx, cy, target, sc, inst, et);
% Project origin to focal plane
origin = topo2inst(origin_topo, cx, cy, target, sc, inst, et);

% Get sweeping directions (if necessary)
if isempty(sweepDir1)
    % Closest polygon side to the spacecraft's ground track position (this
    % will determine the coverage path)
    gt1 = groundtrack(sc, et, target);
    gt2 = groundtrack(sc, et + 500, target);
    gt1 = topo2inst(gt1, cx, cy, target, sc, inst, et);
    gt2 = topo2inst(gt2, cx, cy, target, sc, inst, et + 500);
    [sweepDir1, sweepDir2] = closestSide(gt1, gt2, targetArea, angle);
end

if isempty(topo_grid) % First iteration
    % Focal plane grid discretization
    grid = grid2D(fpref, olapx, olapy, origin, targetArea);
else
    % Replanning Sidewinder: unlike Sidewinder, which only calls
    % planSidewinderTour once, replanning Sidewinder "replans" after
    % every iteration (i.e., footprint). Therefore, in order not to
    % lose the original coverage path direction, a boolean variable is
    % defined to identify when the tour is hopping to a new row/column.
    % If this variable was not put, planSidewinderTour would always
    % start at the right/left or top/bottom side of the grid (depending
    % on the spacecraft's moving direction). So, for example, if we
    % have a horizontal sweep, instead of doing
    % this:
    % ---->-----|
    %           |
    % |---<-----|
    % |
    % |--->------
    % It would be doing this (discontinuities between lines)
    % ----------->
    % >-----------
    % ------------>
    %
    % Likewise, if we have a vertical sweep, instead of doing this:
    % |  ----
    % |  |  |
    % ⌄  ^  ⌄
    % |  |  |
    % |  |  |
    % ----  |
    % It would be doing this (discontinuities between lines)
    % | | |
    % | | |
    % ⌄ ⌄ ⌄
    % | | |
    % | | |
    % | | |
    %
    % The following algorithm identifies if the next element in the
    % tour will imply a change of row/column in the grid
    [i, j] = find(~cellfun('isempty', topo_grid));
    old = []; curr = [];
    for k=1:numel(i)
        if isequal(old_origin_topo, topo_grid{i(k), j(k)})
            old   = [i(k), j(k)];
        elseif isequal(origin_topo, topo_grid{i(k), j(k)})
            curr  = [i(k), j(k)];
        end
        if ~isempty(old) && ~isempty(curr)
            break;
        end
    end
    
    % In case the seed ends up being sterile, we need to retrieve the old 
    % sweeping directions
    oldSweepDir1 = sweepDir1; oldSweepDir2 = sweepDir2;

    % Check for sweepDir changes in the Boustrophedon traversal
    if isequal(sweepDir1, 'east') || isequal(sweepDir1, 'west') % vertical sweep
        if curr(2) ~= old(2)
            % direction must be the contrary to the current one
            if isequal(sweepDir2, 'north'), sweepDir2 = 'south';
            else, sweepDir2 = 'north'; end
        end
    elseif isequal(sweepDir1, 'north') || isequal(sweepDir1, 'south')
        if curr(1) ~= old(1) % change of row between replans
            % next row's coverage path's
            if isequal(sweepDir2, 'west'), sweepDir2 = 'east';
            else, sweepDir2 = 'west'; end
        end
    end

    % Replanning sidewinder: optimize grid origin in order to avoid
    % potential taboo tiles
    grid = optimizeGridOrigin(origin, fpref, olapx, olapy, ...
          targetArea, sweepDir1, sweepDir2);
end

if ~isempty(grid)
    % Boustrophedon decomposition
    itour = boustrophedon(grid, sweepDir1, sweepDir2);

    % Transform coordinates
    topo_grid = inst2topo(grid, cx, cy, target, sc, inst, et);
    topo_tour = inst2topo(itour, cx, cy, target, sc, inst, et);
    % Remove empty elements from the tour, which may result from unobservable
    % regions within the planned plath
    emptyCells = cellfun(@isempty, topo_tour); % find indices of empty cells
    topo_tour(emptyCells) = []; % remove empty cells
    
    % Flag indicating if the seed was sterile
    emptySeed = false;
else
    % Sterile seed (the observation point is no longer worth it)
    topo_tour(1) = [];
    sweepDir1 = oldSweepDir1; sweepDir2 = oldSweepDir2;

    % Flag indicating that current seed was sterile (because it is no
    % longer visible or its allocated tile does not comply with the minimum
    % threshold of ROI coverage)
    emptySeed = true;
end
end