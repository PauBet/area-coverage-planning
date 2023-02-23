function [tour, grid] = planSidewinderTour(target, sc, t, roi, fprint0, ...
    gamma, olapx, olapy, cside)
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
%   > gamma:        origin of the grid, in latitudinal coordinates, in deg
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage
% 
% Outputs:
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables
tour = {}; % list of planned observations
[~, targetFrame, ~] = cspice_cnmfrm(target); % target body-fixed reference 
% frame

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

% 2D grid discretization
grid = grid2D(fprint0, olapx, olapy, gamma, roi);

% The origin of the coverage path depends on the spacecraft ground track
% position
switch cside
    case {'up','down'} % Horizontal sweep
        if counter == 1
            if (sclon - sclon_) >= 0 % sc is moving to the left
                rightsweep = false;
            else % right
                rightsweep = true;
            end
        end

        if rightsweep
            bearing = true; % if the spacecraft is moving right (in the
            % topography map) then the coverage path should start at the
            % position furthest to the left (left -> right direction)
        else
            bearing = false; % if the spacecraft is moving left (in the
            % topography map) then the coverage path should start at the
            % position furthest to the right (right -> left direction)
        end

        for i=1:size(grid,1)
            % Sweep across latitude
            if isequal(cside, 'up')
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
                    y = grid{irow, icol}(2);
                    x = grid{irow, icol}(1);
                    tour{end + 1} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each 
            % row sweeping, i.e. left (highest lon) to right (lowest lon) 
            % or vice versa
        end
        
        % Replanning Sidewinder: unlike Sidewinder, which only calls
        % planSidewinderTour once, replanning Sidewinder "replans"
        % after every iteration (i.e., footprint). Therefore, in order not
        % to lose the original coverage path direction, a boolean variable
        % is defined to identify when the tour is hopping to a new row, in
        % which case the tour must follow the same column. If this variable
        % was not put, planSidewinderTour would always start at the
        % right/left side of the grid (depending on the spacecraft's moving
        % direction).
        % So, for example, instead of doing this:
        % ---->-----|
        %           |
        % |---<-----|
        % |
        % |--->------
        % It would be doing this (discontinuities between lines)
        % ----------->
        % >-----------
        % ------------>
        % The following algorithm identifies if the next element in the
        % tour will imply a change of row in the grid
        if counter > 1 && length(tour) > 1
            flag = false;
            for i=1:size(grid, 1)
                for j=1:size(grid, 2)
                    if ~isempty(grid{i, j})
                        if tour{1}(1) == grid{i,j}(1) && ...
                                tour{1}(2) == grid{i,j}(2)
                            curr = [i, j]; % current depicted element
                        elseif tour{2}(1) == grid{i,j}(1) && ...
                                tour{2}(2) == grid{i,j}(2)
                            new  = [i, j]; % next depicted element
                        end
                    end
                    if exist('new', 'var') && exist('curr', 'var')
                        flag = true;
                        break;
                    end
                end
                if flag
                    break;
                end
            end

            if new(1) ~= curr(1) % change of row between replans
                rightsweep = not(rightsweep); % next row's coverage path's
                % direction must be the contrary to the current one
            end
        end
    case {'right','left'} % Vertical sweep
        if counter == 1
            if (sclat - sclat_) < 0 % up
                downsweep = false;
            else % down
                downsweep = true; 
            end
        end
        
        if downsweep
            bearing = true; % if the spacecraft is moving down (in the
            % topography map) then the coverage path should start at the
            % position furthest to the top (top -> down direction)
        else
            bearing = false; % if the spacecraft is moving up (in the
            % topography map) then the coverage path should start at the
            % position furthest to the bottom (down -> top direction)
        end

        for i=1:size(grid,2)
            % Sweep across longitude
            if isequal(cside, 'left')
                icol = size(grid, 2) - i + 1;
            else
                icol = i;
            end
            for j=1:size(grid, 1)
                if ~bearing
                    irow = j;
                else
                    irow = size(grid, 1) + 1 - j;
                end
                if ~isempty(grid{irow, icol})
                    y = grid{irow, icol}(2);
                    x = grid{irow, icol}(1);
                    tour{end + 1} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each
            % column sweeping, i.e. up (highest lat) to down (lowest
            % lat) or vice versa
        end
        
        % Replanning Sidewinder: unlike Sidewinder, which only calls
        % planSidewinderTour once, replanning Sidewinder "replans"
        % after every iteration (i.e., footprint). Therefore, in order not
        % to lose the original coverage path direction, a boolean variable
        % is defined to identify when the tour is hopping to a new row, in
        % which case the tour must follow the same column. If this variable
        % was not put, planSidewinderTour would always start at the
        % top/bottom side of the grid (depending on the spacecraft's moving
        % direction).
        % So, for example, instead of doing this:
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
        % The following algorithm identifies if the next element in the
        % tour will imply a change of column in the grid
        if counter > 1 && length(tour) > 1
            flag = false;
            for i=1:size(grid, 1)
                for j=1:size(grid, 2)
                    if ~isempty(grid{i, j})
                        if tour{1}(1) == grid{i,j}(1) && ...
                                tour{1}(2) == grid{i,j}(2)
                            curr = [i, j]; % current depicted element
                        elseif tour{2}(1) == grid{i,j}(1) && ...
                                tour{2}(2) == grid{i,j}(2)
                            new  = [i, j]; % next depicted element
                        end
                    end
                    if exist('new', 'var') && exist('curr', 'var')
                        flag = true;
                        break;
                    end
                end
                if flag
                    break;
                end
            end

            if new(2) ~= curr(2) % change of column between replans
                downsweep = not(downsweep); % next row's coverage path's
                % direction must be the contrary to the current one
            end
        end
end

end