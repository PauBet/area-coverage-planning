function [dir1, dir2, curr, old] = sweepdir2(gamma0, gamma, grid, dir1, ...
    dir2)
% Given a discretization grid, the current and previous observation points
% and the original position of the observer's ground track w.r.t. to the
% area of interest, this function outputs the sweeping direction (which
% alternates every column/row of the grid cell matrix). See
% planSidewinderTour for further details.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        sweep = sweepdir(gamma0, gamma, map, cside, sweep)
%
% Inputs:
%   > gamma0:       previous observation point
%   > gamma:        current observation point
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > cside:        string defining the spacecraft ground track position
%                   with respect to the roi, i.e. 'up', 'down', 'left' or
%                   'right'. See closestSide for further details
%   > sweep:        boolean variable that defines the sweeping direction,
%                   depending on the 'cside' variable. If cside is 'up' or
%                   'down' it is a horizontal sweep, if cside is 'left' or
%                   'right' it is a vertical sweep. In this context, the
%                   sweeping direction defines if the direction is
%                   left to right / up to down or viceversa, respectively.
% 
% Outputs:
%   > sweep:        updated boolean variable that defines the sweeping
%                   direction
%   > curr:         current observation point location indices in 'map'
%   > old:          previous observation point location indices in 'map'

% Pre-allocate variables
curr = []; old = [];

% Previous check...
if isempty(gamma0)
    return;
end

% Find the current observation point location in the grid
% Identify if the current observation point implies a change of row/column
% in the grid (change of sweeping direction)
el = find(~cell2mat(cellfun(@(x)any(isnan(x)), grid, ...
    'UniformOutput', false))); % get non-NaN elements in the map
[i, j] = ind2sub(size(grid), el);
for k=1:numel(i)
    if isequal(gamma0, grid{i(k), j(k)})
        old   = [i(k), j(k)];
    elseif isequal(gamma, grid{i(k), j(k)})
        curr  = [i(k), j(k)];
    end
    if ~isempty(old) && ~isempty(curr)
        break;
    end
end

% Check...
if isempty(curr)
    error("Current origin point not found in the grid...")
end

% dir2 sweeping direction definition
% If there is a change in row/column, direction must be the contrary to the
% current one
if ismember(dir1, {'north', 'south'})
    if curr(1) ~= old(1) % change of row between replans
        if isequal(dir2, 'west'), dir2 = 'east';
        else, dir2 = 'west'; end
    end
else
    if curr(2) ~= old(2) % change of column between replans
        if isequal(dir2, 'north'), dir2 = 'south';
        else, dir2 = 'north'; end
    end
end

end