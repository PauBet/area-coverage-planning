function [map, tour] = grid2Dresize(gamma, map, tour, dirx, diry, ...
    neww, newh)
% Given a uniform grid (map) and a referential point in this grid (gamma), 
% this function scales the former to the required grid spacing.
% It also updates the elements in tour accordingly (see FrontierRepair
% function)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        [map, tour] = grid2Dresize(gamma, map, tour, dirx, 
%                   diry, neww, newh)
%
% Inputs:
%   > gamma:        reference point in the grid. The rest of the grid is
%                   going to be re-scaled with respect to this point.
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > dirx:         grid x-axis direction in the lon-lat map
%   > diry:         grid y-axis direction in the lon-lat map
%   > neww:         new grid spacing in the x-direction
%   > newh:         new grid spacing in the y-direction
% 
% Outputs:
%   > map:          updated, re-scaled cell matrix of grid points.
%   > tour:         updated, re-scaled cell matrix of the successive 
%                   planned observations.
%
% Note: the x-y axes of the map do not necessarily follow the row-column
% axes of the cell matrix. This is because the points in the map follow
% shifted x-y axes according to the reference footprint's angle (the 
% reference footprint sets the grid spacing and orientation). See grid2D
% function for further details

% Pre-allocate variables
Nx = size(map, 1);
Ny = size(map, 2);

% Find the original point (gamma) location in the grid
flag = false;
indrow = []; indcol = [];
for indrow=1:Nx
    for indcol=1:Ny
        if ~isempty(map{indrow, indcol})
            if norm(map{indrow, indcol} - gamma) < 1e-5
                flag = true;
                break;
            end
        end
    end
    if flag, break; end
end
if isempty(indrow) || isempty(indcol)
    error("Gamma not found in the map")
end

% Grid resize
for i=1:Nx
    deltarow = i - indrow;
    for j=1:Ny
        deltacol = j - indcol;
        if ~any(isempty(map{i ,j})) && ~any(isnan(map{i, j}))

            % Find element in tour
            flag = false;
            for tt=1:length(tour)
                if norm(map{i, j} - tour{tt}) < 1e-5
                    flag = true;
                    break;
                end
            end

            if flag
                % Resize element
                map{i, j}(1) =  map{indrow, indcol}(1) + ...
                    deltacol*neww*dirx(1) - deltarow*newh*diry(1);
                map{i, j}(2) =  map{indrow, indcol}(2) + ...
                    deltacol*neww*dirx(2) - deltarow*newh*diry(2);
                tour{tt} = map{i, j};
            end
        end
    end
end
end