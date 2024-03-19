function n = getMapNeighbours(varargin)
% Given a grid of points (cell matrix) and an element, this function 
% outputs the neighbouring points of the current point in the matrix.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% Last Rev.:    06/2023
% 
% Usage:        n = getMapNeighbours(indrow, indcol, map)
%               n = getMapNeighbours(indrow, indcol, map, search)
%
% Inputs:
%   > indrow:       int row index of the matrix element (grid point)
%   > indcol:       int column index of the matrix element (grid point)
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > search:       string that defines if the function shall differentiate
%                   between 'cardinal' and 'diagonal' searches. Otherwise, 
%                   the 8 adjacent points are visited
% 
% Outputs:
%   > n:            cell array with the non-NaN neighbouring points in
%                   the map

% Pre-allocate variables
n = {};
indrow = varargin{1};
indcol = varargin{2};
map    = varargin{3};

% Previous checks...
% Searching element is not in the boundaries
if indrow == 1 || indrow == size(map, 1) || indcol == 1 ...
        || indcol == size(map, 2)
    error("Searching element cannot be in the map boundaries");
end
% Future work: check if first and last rows and columns are NaN

% Search neighbors of the given element in the map
if nargin < 4
    aux_n = cell(8,1);
    if ~any(isnan(map{indrow, indcol})) && ~isempty(map{indrow, indcol})
        aux_n{1} = map{indrow - 1, indcol + 1}; % northeast
        aux_n{2} = map{indrow    , indcol + 1}; % east
        aux_n{3} = map{indrow + 1, indcol + 1}; % southeast
        aux_n{4} = map{indrow - 1, indcol    }; % north
        aux_n{5} = map{indrow + 1, indcol    }; % south
        aux_n{6} = map{indrow - 1, indcol - 1}; % northwest
        aux_n{7} = map{indrow    , indcol - 1}; % west
        aux_n{8} = map{indrow + 1, indcol - 1}; % southwest
    end
else
    % cardinal or diagonal search
    switch varargin{4}
        case 'cardinal'
            aux_n{1} = map{indrow - 1, indcol    }; % north
            aux_n{2} = map{indrow    , indcol + 1}; % east
            aux_n{3} = map{indrow + 1, indcol    }; % south
            aux_n{4} = map{indrow    , indcol - 1}; % west
        case 'diagonal'
            aux_n{1} = map{indrow - 1, indcol - 1}; % northwest
            aux_n{2} = map{indrow - 1, indcol + 1}; % northeast
            aux_n{3} = map{indrow + 1, indcol + 1}; % southeast
            aux_n{4} = map{indrow + 1, indcol - 1}; % southwest
    end
end

% Output neighbours (not empty nor NaN)
for i=1:length(aux_n)
    if ~any(isnan(aux_n{i})) && ~isempty(aux_n{i})
        n{end + 1} = aux_n{i};
    end
end

end