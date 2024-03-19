function grid = map2grid(map)
% This function takes a map, defined as a cell matrix where the outermost 
% rows and columns are filled with NaNs to denote boundaries or areas 
% outside of interest, and converts it into a grid by removing these NaN 
% borders. Within the resulting grid, any cell containing NaN values is 
% emptied, signifying that it does not contain useful data or represents an
% area outside the region of interest.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        grid = map2grid(map)
%
% Input:
%   > map:          cell array representing the map. It is the input grid 
%                   with an added border of NaN values. Inside the grid, 
%                   empty cells are replaced with [NaN NaN]
%
% Output:
%   > grid:         cell array that contains a 2-element vector for grid 
%                   points within the ROI, or is empty for points outside 
%                   the ROI or excluded from coverage.

grid = map(2:end-1, 2:end-1);
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        if any(isnan(grid{i, j}), 'all')
            grid{i, j} = [];
        end
    end
end
end