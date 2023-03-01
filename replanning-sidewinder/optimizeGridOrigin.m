function [grid, gamma] = optimizeGridOrigin(gamma0, fp0, olapx, olapy, ...
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
% Usage:        tour = planSidewinderTour(closestSide, roi, fprint0, gamma)
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
%                   is the spacecraft ground track position with respect to
%                   the edges of the target area. See closestSide function
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
deltax = 0.1*fp0.sizex; % displacement value in the x direction
deltay = 0.1*fp0.sizey; % displacement value in the y direction

% if the grid is not optimal nor the displacement value has not reached its
% maximum yet...
while ~opt && abs(gamma(1) - gamma0(1)) <= fp0.sizex/2 && ...
        abs(gamma(2) - gamma0(2)) <= fp0.sizey/2

    % Discretize the non-covered roi space (flood-fill), seeded with gamma
    gamma = gamma + delta;
    [grid, dirx, diry] = grid2D(fp0, olapx, olapy, gamma, roi);

    % The grid shifting may go in both directions (dirx and diry) and, 
    % in that case, is accumulative
    delta = [0 0];
    
    % Find which position does gamma occupy in this grid
    flag = 0;
    for i=1:size(grid,1)
        for j=1:size(grid,2)
            if ~isempty(grid{i,j})
                if norm(grid{i,j} - gamma') < 1e-3
                    ind = [i j];
                    flag = 1;
                    break;
                end
            end
        end
        if flag
            break;
        end
    end

    if exist('ind', 'var')
        switch cside
            case {'up', 'down'} % horizontal sweep

                if isequal(cside, 'down') % spacecraft is towards roi's
                    % bottom
                    if ind(1) > 1 % if gamma is not in the first row of the
                        % grid, move the grid upwards
                        if sum(cellfun(@any, grid(ind(1) - 1, 1:ind(2))))
                            delta = delta + deltay*diry;
                        end
                    end
                else % spacecraft is towards roi's top
                    if ind(1) < size(grid, 1) % if gamma is not in the
                        % last row of the grid, move the grid downwards
                        if sum(cellfun(@any, grid(ind(1) + 1, 1:ind(2))))
                            delta = delta - deltay*diry;
                        end
                    end
                end

                if dir % tour is moving to the right (left -> right dir.)
                    if ind(2) > 1 % if gamma is not in the first column of
                        % the grid, move the grid leftwards
                        if sum(cellfun(@any, grid(1:ind(1), ind(2) - 1)))
                            delta = delta - deltax*dirx;
                        end
                    end
                else % tour is moving to the left (right -> left dir.)
                    if ind(2) < size(grid, 2) % if gamma is not in the last
                        % column of the grid, move the grid rightwards
                        if sum(cellfun(@any, grid(1:ind(1), ind(2) + 1)))
                            delta = delta + deltax*dirx;
                        end
                    end
                end

            case {'right', 'left'} % vertical sweep

                if isequal(cside, 'right')
                    if ind(2) > 1 % if gamma is not in the first column of
                        % the grid, move the grid leftwards
                        if sum(cellfun(@any, grid(1:ind(1), ind(2) - 1)))
                            delta = delta - deltax*dirx;
                        end
                    end
                else
                    if ind(2) < size(grid, 2) % if gamma is not in the
                        % last column of the grid, move the grid rightwards
                        if sum(cellfun(@any, grid(1:ind(1), ind(2) + 1)))
                            delta = delta + deltax*dirx;
                        end
                    end
                end

                if dir % downsweep = true. tour is moving to the bottom
                    if ind(1) > 1 % if gamma is not in the first row of
                        % the grid, move the grid upwards
                        if sum(cellfun(@any, grid(ind(1) - 1, 1:ind(2))))
                            delta = delta + deltay*diry;
                        end
                    end
                else % tour is moving to the top
                    if ind(1) < size(grid, 1) % if gamma is not in the last
                        % row of the grid, move the grid downwards
                        if sum(cellfun(@any, grid(ind(1) + 1, 1:ind(2))))
                            delta = delta - deltay*diry;
                        end
                    end
                end
        end
    else
        error("gamma (grid origin or seed) not found in the grid");
    end
    
    % If delta = [0 0], no taboo tiles were found
    if norm(delta) == 0
        opt = true;
    end
end
end