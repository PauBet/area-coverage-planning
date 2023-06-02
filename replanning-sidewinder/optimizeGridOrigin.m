function [grid, gamma] = optimizeGridOriginvf(gamma0, fp0, olapx, olapy, ...
    roi, dir, cside)
% This function optimizes the grid origin of an uncovered ROI area. The
% grid origin is the next tile in the tour that is going to be observed.
% Sometimes, within the replanning Sidewinder algorithm, the allocated cell
% for the previous observation was bigger than the actual footprint
% size... therefore, some part of the allocated cell was actually not
% covered and was left for the next observation. If the
% aforementioned uncovered area is significant, the algorithm is going to
% re-program that observation (at the same spot), a.k.a. taboo tile,
% possibly triggering a grid lock. This function intends to prevent that by
% shifting the grid backwards so the next observation may cover (partially
% or totally) the previously uncovered area, but not completely going
% backwards (to the point where the next and previous observations are
% equal)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [grid, gamma] = optimizeGridOriginv2(gamma0, fp0, olapx, 
%    olapy, roi, dir, cside)
%
% Inputs:
%   > gamma0:       origin of the grid, in latitudinal coordinates, in deg
%   > fprint0:      struct containing the footprint parameters that are
%                   going to be used to define the grid discretization
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage
%   > roi:          matrix containing the roi of the ROI polygon. The
%                   vertex points are expressed in 2D. 
%       # roi(:,1) correspond to the x values of the roi
%       # roi(:,2) correspond to the y values of the roi
%   > dir:          boolean variable that states in which direction is the
%                   coverage path currently touring the grid
%   > cside:        given a region-of-interest, this function defines what 
%                   is the spacecraft's ground track position with respect
%                   to the edges of the target area. See closestSide 
%                   function
% 
% Outputs:
%   > grid:         Flood-fill grid discretization of the 
%                   region-of-interest given a reference footprint(fprint0)
%

% Variables
gamma = gamma0; % grid origin (seed), i.e., next planned observation in the
% tour
opt = false; % boolean that states if the grid is optimal -i.e., no taboo 
% tiles- (opt = true) or not (opt = false)
delta = [0 0]; % array that defines the displacement of the grid
deltax = 0.1*fp0.sizex; % displacement value in x direction
deltay = 0.1*fp0.sizey; % displacement value in y direction

% while the grid is not optimal nor the displacement value has reached its
% maximum yet...
grid = [];
while ~opt && abs(gamma(1) - gamma0(1)) <= fp0.sizex/4 && ...
        abs(gamma(2) - gamma0(2)) <= fp0.sizey/4

    % Discretize the non-covered roi space (flood-fill), seeded with gamma
    gamma = gamma + delta;
    [grid, dirx, diry] = grid2D(fp0, olapx, olapy, gamma, roi);

    % In case the original gamma0 was outside the roi's area (and its 
    % allocated cell was inferior to the estabished threshold in order to
    % add it to the grid) -> move gamma towards the centroid to optimize
    % its area coverage
    if isempty(grid)
        [cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
        dirc = [cx - gamma0(1), cy - gamma0(2)];
        dirc = dirc/norm(dirc);
        delta = delta + [deltax*dirc(1), deltay*dirc(2)];
        continue
    end

    % The grid shifting may go in both directions (dirx and diry) and,
    % in that case, is accumulative
    delta = [0 0];
    
    % Find which position does gamma occupy in this grid
    for i=1:numel(grid)
        if ~isempty(grid{i}) && norm(grid{i} - gamma') < 1e-3
            [ind_row, ind_col] = ind2sub(size(grid), i);
            break;
        end
    end

    if exist('ind_row', 'var')
        switch cside
            case {'up', 'down'} % horizontal sweep
                
                % Y direction
                if strcmp(cside, 'down') % spacecraft is towards roi's
                    % bottom
                    if ind_row > 1 && any(cellfun(@any, ...
                            grid(ind_row - 1, :)), 'all')

                        % Move the grid "upwards"
                        delta = delta + deltay*diry;
                    end

                else % spacecraft is towards roi's top
                    if ind_row < size(grid, 1) && any(cellfun(@any, ...
                            grid(ind_row + 1, :)), 'all') 
                        
                        % Move the grid "downwards"
                        delta = delta - deltay*diry;
                    end
                end
                
                % X direction
                if dir % tour is moving to the right (left -> right dir.)
                    if ind_col > 1 && any(cellfun(@any, grid(ind_row, ...
                            ind_col - 1)), 'all')

                        % Move the grid "leftwards"
                        delta = delta - deltax*dirx;
                    end
                else % tour is moving to the left (right -> left dir.)
                    if ind_col < size(grid, 2) && any(cellfun(@any, ...
                            grid(ind_row, ind_col + 1)), 'all')

                        % Move the grid "rightwards"
                        delta = delta + deltax*dirx;
                    end
                end

            case {'right', 'left'} % vertical sweep
                
                % X direction
                if strcmp(cside, 'right')
                    if ind_col > 1 && any(cellfun(@any, ...
                            grid(:, ind_col - 1)), 'all')

                        % Move the grid leftwards
                        delta = delta - deltax*dirx;
                    end
                else
                    if ind_col < size(grid, 2) && any(cellfun(@any, ...
                            grid(:, ind_col + 1)), 'all')
                        
                        % Move the grid rightwards
                        delta = delta + deltax*dirx;
                    end
                end
                
                % Y direction
                if dir % downsweep = true. tour is moving to the bottom
                    if ind_row > 1 && any(cellfun(@any, ...
                            grid(ind_row - 1, ind_col)), 'all')

                        % Move the grid upwards
                        delta = delta + deltay*diry;
                    end
                else % tour is moving to the top
                    if ind_row < size(grid, 1) && any(cellfun(@any, ...
                            grid(ind_row + 1, ind_col)), 'all')

                        % Move the grid downwards
                        delta = delta - deltay*diry;
                    end
                end
        end
    end

    % If delta = [0 0], no taboo tiles (worth correcting) were found
    if norm(delta) == 0
        opt = true;
        gamma0 = gamma;
    end
end

if isempty(grid)
    return;
end

% Check if there are any taboo tiles that should be deleted... (at the
% expense of potential uncovered area)
% This matrix points at the taboo tiles of grid (taboo_idx(i, j) = 1 if
% there is a taboo tile at element [i, j] of the grid)
taboo_idx = zeros(size(grid));

switch cside
    case {'up', 'down'} % horizontal sweep

        if strcmp(cside, 'down') % spacecraft is towards roi's
            % bottom
            if ind_row > 1 && any(cellfun(@any, ...
                    grid(1:(ind_row - 1), :)), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid(1:(ind_row - 1), :)));
                taboo_idx(row, col) = 1;
            end
        else % spacecraft is towards roi's top
            if ind_row < size(grid, 1) && any(cellfun(@any, ...
                    grid((ind_row + 1):end, :)), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid((ind_row + 1):end, :)));
                row = row + ind_row; % adjust row indices to be
                % relative to the entire grid
                taboo_idx(row, col) = 1;
            end
        end

        if dir % tour is moving to the right (left -> right dir.)
            if ind_col > 1
                if isequal(cside, 'down') && any(cellfun(@any, ...
                        grid(1:ind_row, 1:(ind_col - 1))), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:ind_row, 1:(ind_col - 1))));
                    taboo_idx(row, col) = 1;
                elseif isequal(cside, 'up') && any(cellfun(@any, ...
                        grid(ind_row:end, 1:(ind_col - 1))), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(ind_row:end, 1:(ind_col - 1))));
                    row = row + ind_row - 1; % adjust row indices to be
                    % relative to the entire grid
                    taboo_idx(row, col) = 1;
                end
            end
        else % tour is moving to the left (right -> left dir.)
            if ind_col < size(grid, 2)
                if isequal(cside, 'down') && any(cellfun(@any, ...
                        grid(1:ind_row, (ind_col + 1):end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:ind_row, (ind_col + 1):end)));
                    col = col + ind_col; % adjust column indices to
                    % be relative to the entire grid
                    taboo_idx(row, col) = 1;
                elseif isequal(cside, 'up') && any(cellfun(@any, ...
                        grid(ind_row:end, (ind_col + 1):end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(ind_row:end, (ind_col + 1):end)));
                    col = col + ind_col; % adjust column indices to
                    % be relative to the entire grid
                    row = row + ind_row - 1; % adjust row indices
                    % to be relative to the entire grid
                    taboo_idx(row, col) = 1;
                end

            end
        end

    case {'right', 'left'} % vertical sweep

        if strcmp(cside, 'right')
            if ind_col > 1 && any(cellfun(@any, ...
                    grid(:, 1:(ind_col - 1))), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid(:, 1:(ind_col - 1))));
                taboo_idx(row, col) = 1;
            end
        else
            if ind_col < size(grid, 2) && any(cellfun(@any, ...
                    grid(:, (ind_col + 1):end)), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid(:, (ind_col + 1):end)));
                col = col + ind_col; % adjust column indices to
                % be relative to the entire grid
                taboo_idx(row, col) = 1;
            end
        end

        if dir % downsweep = true. tour is moving to the bottom
            if ind_row > 1
                if isequal(cside, 'right') && ...
                        any(cellfun(@any, grid(1:(ind_row - 1), ...
                        1:ind_col)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:(ind_row - 1), 1:ind_col)));
                    taboo_idx(row, col) = 1;

                elseif isequal(cside, 'left') && ...
                        any(cellfun(@any, grid(1:(ind_row - 1), ...
                        ind_col:end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:(ind_row - 1), ind_col:end)));
                    col = col + ind_col - 1;
                    taboo_idx(row, col) = 1;
                end
            end
        else % tour is moving to the top
            if ind_row < size(grid, 1)

                if isequal(cside, 'right') && ...
                        any(cellfun(@any, grid((ind_row + 1):end, ...
                        1:ind_col)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid((ind_row + 1):end, 1:ind_col)));
                    row = ind_row + row;
                    taboo_idx(row, col) = 1;
                elseif isequal(cside, 'left') && ...
                        any(cellfun(@any, grid((ind_row + 1):end, ...
                        ind_col:end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid((ind_row + 1):end, ind_col:end)));
                    row = ind_row + row;
                    col = ind_col + col - 1;
                    taboo_idx(row, col) = 1;
                end
            end
        end
end

% Delete taboo tiles
grid(taboo_idx == 1) = {[]};

% After deleting the taboo tiles, check if there are empty rows and/or
% columns and erase those from the grid...
empty_rows    = all(cellfun('isempty', grid), 2);
empty_columns = all(cellfun('isempty', grid), 1);
grid(empty_rows, :) = []; grid(:, empty_columns) = [];
end