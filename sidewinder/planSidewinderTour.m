function [tour, grid] = planSidewinderTour(target, sc, t, roi, fprint0, ...
    olapx, olapy, cside, tour0, grid0)
% This function discretizes and computes the coverage path of a certain
% ROI in order to build a mosaic image, adapted from [1]. Note that the
% coverage path is oriented as the width direction of the footprint
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        tour = planSidewinderTour(closestSide, roi, fprint0, gamma)
%
% Inputs:
%   > target:       SPICE string name of the target body
%   > sc:           SPICE string name of the spacecraft
%   > t:            time in TDB seconds past J2000 epoch
%   > roi:          matrix containing the roi of the ROI polygon. The
%                   vertex points are expressed in 2D. 
%       # roi(:,1) correspond to the x values of the roi
%       # roi(:,2) correspond to the y values of the roi
%   > fprint0:      struct containing the footprint parameters that are
%                   going to be used to define the grid discretization
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage
%   > cside:        given a region-of-interest, this function defines what 
%                   is the spacecraft ground track position with respect to
%                   the edges of the target area. See closestSide function
%   > tour0:        seeding tour (set of previously planned observations)
%   > grid0:        seeding grid (previously discretized roi)
% 
% Outputs:
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > grid:         cell matrix that contains the grid discretization of 
%                   the region-of-interest
%
% Note: if it is the initial planning -i.e., there is no former tour- 
% grid0 can be empty, and tour may contain just one point (seed) inside 
% the roi (for the grid discretization)
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Previous checks...
if isempty(tour0)
    error("Seeding tour ('tour0') cannot be empty. Input at least " + ...
        "a seeding point for the grid discretization")
end

% Pre-allocate variables
[~, targetFrame, ~] = cspice_cnmfrm(target); % target body-fixed reference 
% frame
curr = []; old = [];
grid = {}; tour = {};

% ReplanningSidewinder variables:
persistent counter; % variable that counts the number of times this
% function is called
persistent downsweep; % boolean that defines the vertical sweep direction 
% (and persists between replans). true if the coverage path is going down,
% otherwise is going up
persistent rightsweep; % boolean that defines the horizontal sweep 
% direction (and persists between replans). true if coverage is going
% right, otherwise is going left

% count number of times this function has been called
if isempty(counter)
    counter = 1;
else
    counter = counter + 1;
end

% Calculate how is the spacecraft ground track position moving along the
% map. The coverage path does not only depend on the spacecraft position
% itself but also its velocity direction
subobs = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t,...
    targetFrame, 'NONE', sc);
[~, sclon, sclat] = cspice_reclat(subobs); % latitudinal coordinates
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr; % [rad] to [deg]
subobs_ = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t + 60,...
    targetFrame, 'NONE', sc);
[~, sclon_, sclat_] = cspice_reclat(subobs_); % latitudinal coordinates
sclon_ = sclon_*cspice_dpr; sclat_ = sclat_*cspice_dpr; % [rad] to [deg]

% 2D grid discretization and setting coverage path's sense of sweeping
while isempty(grid) && length(tour0) > 1
    gamma_old = tour0{1};
    tour0(1) = [];
    gamma    = tour0{1};

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
    [i, j] = find(~cellfun('isempty', grid0));
    for k=1:numel(i)
        if isequal(gamma_old', grid0{i(k), j(k)})
            old   = [i(k), j(k)];
        elseif isequal(gamma', grid0{i(k), j(k)})
            curr  = [i(k), j(k)];
        end
        if ~isempty(old) && ~isempty(curr)
            break;
        end
    end

    if ismember(cside, {'up', 'down'})
        if curr(1) ~= old(1) % change of row between replans
            rightsweep = not(rightsweep); % next row's coverage path's
            % direction must be the contrary to the current one
        end

        % Replanning sidewinder: optimize grid origin in order to avoid
        % potential taboo tiles
        grid = optimizeGridOrigin(gamma, fprint0, olapx, olapy, ...
            roi, rightsweep, cside);
    else
        if curr(2) ~= old(2) % change of row between replans
            downsweep = not(downsweep); % next row's coverage path's
            % direction must be the contrary to the current one
        end

        % Replanning sidewinder: optimize grid origin in order to avoid
        % potential taboo tiles
        grid = optimizeGridOrigin(gamma, fprint0, olapx, olapy, ...
            roi, downsweep, cside);
    end

    curr = []; old = [];
end

% If grid is still empty: end observations
if isempty(grid)
    return;
end

%% Plan tour over the grid discretization
% The origin of the coverage path depends on the spacecraft ground track
% position
switch cside
    case {'up','down'} % Horizontal sweep

        if counter == 1
            if (sclon - sclon_) >= 0 % sc is moving leftwards
                rightsweep = false; % if the spacecraft is moving left (in 
                % the topography map) then the coverage path should start 
                % at the position furthest to the right (right -> left 
                % direction)
            else % sc is moving to the right
                rightsweep = true; % if the spacecraft is moving right (in 
                % the topography map) then the coverage path should start 
                % at the position furthest to the left (left -> right 
                % direction)
            end
        end

        if rightsweep
            bearing = true; % left -> right
        else
            bearing = false; % right -> left
        end
        
        tour = cell(1, nnz(~cellfun('isempty', grid))); % list of planned 
        % observations
        ii = 0;
        for i=1:size(grid,1)
            % Sweep across latitude
            if isequal(cside, 'down')
                irow = i;
            else
                irow = size(grid, 1) - i + 1;
            end
            for j=1:size(grid, 2)
                if ~bearing
                    icol = size(grid, 2) - j + 1;
                else
                    icol = j;
                end
                if ~isempty(grid{irow, icol})
                    ii = ii + 1;
                    y = grid{irow, icol}(2);
                    x = grid{irow, icol}(1);
                    tour{ii} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each 
            % row sweeping, i.e. left (highest lon) to right (lowest lon) 
            % or vice versa
        end
 
    case {'right','left'} % Vertical sweep

        if counter == 1
            if (sclat - sclat_) < 0 % sc is moving upwards
                downsweep = false; % if the spacecraft is moving up (in the
                % topography map) then the coverage path should start at 
                % the position furthest to the bottom (down -> top 
                % direction)
            else % sc is moving down
                downsweep = true;% if the spacecraft is moving down (in the
                % topography map) then the coverage path should start at 
                % the position furthest to the top (top -> down direction)
            end
        end
        
        if downsweep
            bearing = true; % top -> down
        else
            bearing = false; % down -> top
        end
        
        tour = cell(1, nnz(~cellfun('isempty', grid))); % list of planned 
        % observations
        ii = 0;
        for i=1:size(grid,2)
            % Sweep across longitude
            if isequal(cside, 'left')
                icol = size(grid, 2) - i + 1;
            else
                icol = i;
            end
            for j=1:size(grid, 1)
                if bearing
                    irow = j;
                else
                    irow = size(grid, 1) + 1 - j;
                end
                if ~isempty(grid{irow, icol})
                    ii = ii + 1;
                    y = grid{irow, icol}(2);
                    x = grid{irow, icol}(1);
                    tour{ii} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each
            % column sweeping, i.e. up (highest lat) to down (lowest
            % lat) or vice versa
        end
end

end