function [sweep, curr, old] = sweepdir(gamma0, gamma, map, cside, sweep)
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

% Previous checks...
if isempty(sweep)
    error("Sweeping direction is a boolean variable, it cannot be empty")
end

% Find the current observation point location in the grid
% Identify if the current observation point implies a change of row/column
% in the grid (change of sweeping direction)
el = find(~cell2mat(cellfun(@(x)any(isnan(x)), map, ...
    'UniformOutput', false))); % get non-NaN elements in the map
[i, j] = ind2sub(size(map), el);
curr = []; old = [];
for k=1:numel(i)
    if norm(gamma - map{i(k), j(k)}) < 1e-5
        curr = [i(k), j(k)];
    end
    if ~isempty(curr)
        break;
    end
end
if isempty(curr)
    error("Current observation point not found in the map")
end

% If there was no previous observation (empty), no need to check
% the sweeping direction
if isempty(gamma0)
    return;
else
    for k=1:numel(i)
        if norm(gamma0 - map{i(k), j(k)}) < 1e-5
            old = [i(k), j(k)];
        end
        if ~isempty(old)
            break;
        end
    end
end

if ismember(cside, {'up', 'down'})
    if curr(1) ~= old(1) % change of row between replans
        sweep = not(sweep);
    end
else
    if curr(2) ~= old(2) % change of row between replans
        sweep = not(sweep);
    end
end

end