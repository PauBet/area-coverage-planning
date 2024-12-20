function [N, Nind] = checkTaboo(N, Nind, map, ind_row, ind_col, indir1, indir2)
% This function evaluates each tile in a list of potential observation 
% points (N, Nind) and determines whether it should be considered taboo, 
% i.e., unsuitable for selection based on its location relative to a specified
% direction of movement (dir1, dir2) and starting point (ind_row, ind_col).
% Taboo tiles are those that do not conform to the expected movement 
% pattern across the grid
% [Future work]: get rid of boustrophedonMod function
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [N, Nind] = checkTaboo2(N, Nind, map, ind_row, ind_col, dir1, dir2)
%
% Inputs:
%   > N:            cell array containing the values of potential 
%                   observation points
%   > Nind:         cell array containing the indices of potential 
%                   observation points within the map
%   > map:          cell matrix of grid points, where first and last rows 
%                   and columns are NaN to denote boundaries
%   > ind_row:      starting row index for evaluating taboo conditions
%   > ind_col:      starting column index for evaluating taboo conditions
%   > dir1:         primary movement direction in the path plan
%   > dir2:         secondary movement direction, orthogonal to dir1
%
% Outputs:
%   > N, Nind:      updated lists
% 

persistent pdir1;
persistent pdir2;
if isempty(pdir1)
    pdir1 = indir1;
    pdir2 = indir2;
end
% % Previous checks...
% if isempty(N) || isempty(Nind)
%     return;
% end

% Pre-allocate variables
nindel = [];
grid = map2grid(map);
[dir1, dir2] = boustrophedonMod(grid, indir1, indir2);
dirChange = false;
if ~isequal(pdir2, indir2)
    dirChange = true;
end

% Previous checks...
if isempty(N) || isempty(Nind)
    return;
end

% Define the grid boundaries (ending rows and columns in the map)
for i=size(map, 1):-1:1
    el = find(~cell2mat(cellfun(@(x)any(isnan(x)), map(i, :), ...
    'UniformOutput', false)), 1); % get non-NaN elements in the map
    if ~isempty(el)
        Nrow = i;
        break;
    end
end
for i=size(map, 2):-1:1
    el = find(~cell2mat(cellfun(@(x)any(isnan(x)), map(:, i), ...
    'UniformOutput', false)), 1); % get non-NaN elements in the map
    if ~isempty(el)
        Ncol = i;
        break;
    end
end

% Define the grid boundaries (starting rows and columns in the map)
for i=1:size(map, 1)
    el = find(~cell2mat(cellfun(@(x)any(isnan(x)), map(i, :), ...
    'UniformOutput', false)), 1); % get non-NaN elements in the map
    if ~isempty(el)
        Orow = i;
        break;
    end
end
for i=1:size(map, 2)
    el = find(~cell2mat(cellfun(@(x)any(isnan(x)), map(:, i), ...
    'UniformOutput', false)), 1); % get non-NaN elements in the map
    if ~isempty(el)
        Ocol = i;
        break;
    end
end

% Taboo tiles search
for i=1:length(Nind)
    taboo = false;
    indel = Nind{i};

    switch dir1
        case {'north', 'south'} % horizontal sweep

            if strcmp(dir1, 'south') % spacecraft is towards roi's
                % bottom
                if indel(1) < ind_row
                    taboo = true;
                end
            else % spacecraft is towards roi's top
                if indel(1) > ind_row
                    taboo = true;
                end
            end

            if isequal(dir2, 'west') % tour is moving to the right (left -> right dir.)
                if indel(1) == ind_row && indel(2) > ind_col
                    if dirChange && ind_col < Ncol
                        taboo = true;
                    elseif ~dirChange
                        taboo = true;
                    end
                end
            else % tour is moving to the left (right -> left dir.)
                if indel(1) == ind_row && indel(2) < ind_col 
                    if dirChange && ind_col > Ocol
                        taboo = true;
                    elseif ~dirChange
                        taboo = true;
                    end
                end
            end

        case {'east', 'west'} % vertical sweep

            if strcmp(dir1, 'east')
                if indel(2) < ind_col
                    taboo = true;
                end
            else
                if indel(2) > ind_col
                    taboo = true;
                end
            end

            if isequal(dir2, 'south') % downsweep = true. tour is moving to the bottom
                if indel(2) == ind_col && indel(1) < ind_row 
                    if dirChange && ind_row > Orow
                        taboo = true;
                    elseif ~dirChange
                        taboo = true;
                    end
                end
            else % tour is moving to the top
                if indel(2) == ind_col && indel(1) > ind_row
                    if dirChange && ind_row < Nrow
                        taboo = true;
                    elseif ~dirChange 
                        taboo = true;
                    end
                end
            end
    end

    if taboo
        nindel = [nindel i];
    end
end

% Delete taboo tiles
N(nindel) = [];
Nind(nindel) = [];
end

function [currdir1, currdir2] = boustrophedonMod(grid, dir1, dir2)

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
if isequal(dir2, 'east') || isequal(dir2, 'south')
    sweep = true;
else
    sweep = false;
end
currdir1 = dir1;
currdir2 = dir2;

% Plan tour over the grid discretization
% The origin of the coverage path depends on the spacecraft ground track
% position
start = false;
switch dir1
    case {'north','south'} % Horizontal sweep
        
        if sweep
            bearing = true; % left -> right
        else
            bearing = false; % right -> left
        end

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
                    start = true;
                    break;
                end
            end
            if start, break; end
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
                    start = true;
                    break;
                end
            end
            if start, break; end
            bearing = not(bearing); % Switch coverage direction after each
            % column sweeping, i.e. up (highest lat) to down (lowest
            % lat) or vice versa
        end
end

if bearing 
    if isequal(dir2, 'west'), currdir2 = 'east';
    elseif isequal(dir2, 'north'), currdir2 = 'south';
    end
else
    if isequal(dir2, 'east'), currdir2 = 'west';
    elseif isequal(dir2, 'south'), currdir2 = 'north';
    end
end
end

% [Future work]: get rid of boustrophedonMod
% % Shift directions
% if isequal(dir2, 'west'), dir2 = 'east';
% elseif isequal(dir2, 'east'), dir2 = 'west';
% elseif isequal(dir2, 'north'), dir2 = 'south';
% elseif isequal(dir2, 'south'), dir2 = 'north'; 
% end
% 
% % Identify in which direction are we moving
% if ismember(dir1, {'north', 'south'})
%     if isequal(dir2, 'west')
%         if mod(size(map, 1), 2)
%             if ~mod(ind_row-1, 2), dir2 = 'east'; end
%         else
%             if mod(ind_row-1, 2), dir2 = 'east'; end
%         end
%     else
%         if mod(size(map, 1), 2)
%             if ~mod(ind_row-1, 2), dir2 = 'west'; end
%         else
%             if mod(ind_row-1, 2), dir2 = 'west'; end
%         end
%     end
% else
%     if isequal(dir2, 'north')
%         if mod(size(map, 2), 2)
%             if ~mod(ind_col-1, 2), dir2 = 'south'; end
%         else
%             if mod(ind_col-1, 2), dir2 = 'south'; end
%         end
%     else
%         if mod(size(map, 2), 2)
%             if ~mod(ind_col-1, 2), dir2 = 'north'; end
%         else
%             if mod(ind_col-1, 2), dir2 = 'north'; end
%         end
%     end
% end