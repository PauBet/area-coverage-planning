function map = grid2map(grid)
% This function creates a map from a given grid. It adds a border of NaN 
% values around the entire grid. Within the grid, any cells that are empty 
% or have been excluded from observation (because their coverage is deemed 
% too small) are also filled with NaN values
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        map = grid2map(grid)
%
% Input:
%   > grid:         cell array that contains a 2-element vector for grid 
%                   points within the ROI, or is empty for points outside 
%                   the ROI or excluded from coverage.
%
% Output:
%   > map:          cell array representing the map. It is the input grid 
%                   with an added border of NaN values. Inside the grid, 
%                   empty cells are replaced with [NaN NaN]

% Initialize map with additional rows and columns to place the NaN
% boundaries
map = cell(size(grid,1) + 2, size(grid,2) + 2);

% Populate map with data from the grid and NaN borders
for i=1:size(map,1)
    for j=1:size(map,2)
        if i == 1 || j == 1 || i == size(map,1) || ...
                j == size(map,2) || any(isempty(grid{i-1, j-1}))
            map{i,j} = [NaN NaN];
        else
            map{i, j} = grid{i-1, j-1};
        end
    end
end
end