function [N, Nind] = checkTaboo(N, Nind, gamma, map, cside, dir)
% Given a grid discretization, an observation point and a tour direction of
% the grid (coverage path), this function determines if a set of potential
% points (not included in the grid) would potentially become taboo tiles
% (i.e., the inclusion of these observation points could imply going
% backwards in the coverage path)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        [N, Nind] = checkTaboo(N, Nind, gamma, map, dir, cside)
%
% Inputs:
%   > N:            set of potential observation points
%   > Nind:         index matrices of 'N' in the grid ('map')
%   > gamma:        reference point in the grid (next observation)
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > cside:        string defining the spacecraft ground track position
%                   with respect to the roi, i.e. 'up', 'down', 'left' or
%                   'right'. See closestSide for further details
%   > dir:          boolean variable that defines the sweeping direction,
%                   depending on the 'cside' variable. If cside is 'up' or
%                   'down' it is a horizontal sweep, if cside is 'left' or
%                   'right' it is a vertical sweep. In this context, the
%                   sweeping direction defines if the direction is
%                   left to right / up to down or viceversa, respectively.
% 
% Outputs:
%   > N:            set of updated potential observation points, with no
%                   taboo tiles
%   > Nind:         index matrices of 'N' in the grid, with no taboo tiles

% Previous checks...
if isempty(N) || isempty(Nind)
    return;
end

% Pre-allocate variables
ind_row = []; ind_col = [];
for i=1:numel(map)
    if norm(map{i} - gamma) < 1e-5
        [ind_row, ind_col] = ind2sub(size(map), i);
    end
end
if isempty(ind_row) || isempty(ind_col)
    error("Observation point not found in the map")
end
nindel = [];

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

% Taboo tiles search
for i=1:length(Nind)
    taboo = false;
    indel = Nind{i};

    switch cside
        case {'up', 'down'} % horizontal sweep

            if strcmp(cside, 'down') % spacecraft is towards roi's
                % bottom
                if ind_row > 1 && indel(1) < ind_row
                    taboo = true;
                end
            else % spacecraft is towards roi's top
                if ind_row < Nrow && indel(1) > ind_row
                    taboo = true;
                end
            end

            if dir % tour is moving to the right (left -> right dir.)
                if ind_col > 1
                    if isequal(cside, 'down') && indel(1) <= ind_row ...
                            && indel(2) < ind_col
                        taboo = true;
                    elseif isequal(cside, 'up') && indel(1) >= ind_row ...
                            && indel(2) < ind_col
                        taboo = true;
                    end
                end
            else % tour is moving to the left (right -> left dir.)
                if ind_col < Ncol
                    if isequal(cside, 'down') && indel(1) <= ind_row && ...
                            indel(2) > ind_col
                        taboo = true;
                    elseif isequal(cside, 'up') && indel(1) >= ind_row ...
                            && indel(2) > ind_col
                        taboo = true;
                    end

                end
            end

        case {'right', 'left'} % vertical sweep

            if strcmp(cside, 'right')
                if ind_col > 1 && indel(2) < ind_col
                    taboo = true;
                end
            else
                if ind_col < Ncol && indel(2) > ind_col
                    taboo = true;
                end
            end

            if dir % downsweep = true. tour is moving to the bottom
                if ind_row > 1
                    if isequal(cside, 'right') && indel(1) < ind_row && ...
                            indel(2) <= ind_col
                        taboo = true;
                    elseif isequal(cside, 'left') && indel(1) < ind_row ...
                            && indel(2) >= ind_col
                        taboo = true;
                    end
                end
            else % tour is moving to the top
                if ind_row < Nrow
                    if isequal(cside, 'right') && indel(1) > ind_row && ...
                            indel(2) <= ind_col
                        taboo = true;
                    elseif isequal(cside, 'left') && indel(1) > ind_row ...
                            && indel(2) >= ind_col
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