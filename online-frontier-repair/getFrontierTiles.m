function [frontier, indel] = getFrontierTiles(map)
% Given a grid of points (cell matrix) and an element, this function 
% outputs the set of points that have less than 8 neighbours in the grid
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [frontier, indel] = getFrontierTiles(map, tour)
%
% Inputs:
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > gamma0:       previous observation point. This is the exception
%                   element, i.e., 
% 
% Outputs:
%   > frontier:     cell array that contains the frontier tiles in the map
%   > indel:        cell array that contains the indices where the 
%                   frontier tiles are located in 'map'

% Pre-allocate variables
frontier = {};
indel    = {};

for i=1:numel(map)
    % For each observation point in the discretized grid...
    [j, k] = ind2sub(size(map), i);
    tile   = map{j, k};

    if ~any(isnan(tile))
        % Get the observation neighbours of the current observation points
        n = getMapNeighbours(j, k, map);

        % If the observation point has less than 8 planned neighbours, then 
        % it is regarded as a frontier tile
        if length(n) < 8
            frontier{end+1} = tile;
            indel{end+1}    = [j, k];
        end
    end
end
end