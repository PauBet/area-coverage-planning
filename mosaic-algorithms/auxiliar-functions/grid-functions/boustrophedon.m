function tour = boustrophedon(grid, dir1, dir2)
% This function plans an observation tour over a specified grid, creating a
% path that covers the area in alternating rows/columns, according to a 
% specified input direction
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        tour = boustrophedon(grid, dir1, dir2)
%
% Inputs:
%   > grid:        cell array where each cell contains the coordinates of
%                  an observation point or is empty if there is no point
%   > dir1:        primary direction of the sweep
%   > dir2:        secondary direction of the sweep. This defines if the
%                  traversal is going to be performed either in 
%                  alternating rows or columns
% 
% Outputs:
%   > tour:        ordered cell array of points representing the planned 
%                  tour. Each element of the array is a 2-element vector 
%                   indicating a point on the grid to be observed

% Previous check...
if isequal(dir1, 'north') || isequal(dir1, 'south')
    if ~isequal(dir2, 'east') && ~isequal(dir2, 'west')
        error("Sweeping direction is not well defined")
    end
elseif isequal(dir1, 'east') || isequal(dir1, 'west')
    if ~isequal(dir2, 'north') && ~isequal(dir2, 'south')
        error("Sweeping direction is not well defined")
    end
end

% Pre-allocate variables
tour = {};
% Particular case: length(grid) = 1
if length(grid) == 1 
    tour = grid(1);
    return;
end
if isequal(dir2, 'east') || isequal(dir2, 'south')
    sweep = true;
else
    sweep = false;
end

%% Plan tour over the grid discretization
% The origin of the coverage path depends on the spacecraft ground track
% position
switch dir1
    case {'north','south'} % Horizontal sweep

        if sweep
            bearing = true; % left -> right
        else
            bearing = false; % right -> left
        end

        tour = cell(1, nnz(~cellfun('isempty', grid))); % list of planned
        % observations
        ii = 0;
        for i=1:size(grid,1)
            % Sweep across latitude
            if isequal(dir1, 'south')
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

    case {'east','west'} % Vertical sweep

        if sweep
            bearing = true; % top -> down
        else
            bearing = false; % down -> top
        end

        tour = cell(1, nnz(~cellfun('isempty', grid))); % list of planned
        % observations
        ii = 0;
        for i=1:size(grid,2)
            % Sweep across longitude
            if isequal(dir1, 'west')
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