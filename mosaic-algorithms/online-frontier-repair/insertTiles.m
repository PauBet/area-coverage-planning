function map = insertTiles(varargin)
% This function includes new observation points in a planned tour
% observations
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        [tour, map] = insertTiles(map, newp, indp)
%
% Inputs:
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > newp:         cell matrix of the new observation points to be
%                   included in 'tour'
%   > indp:         cell matrix of the new observation points locations to
%                   be included in 'map'
% 
% Outputs:
%   > map:         updated cell matrix of grid points

% Pre-allocate variables
map  = varargin{1};
newp = varargin{2};
indp = varargin{3};

%% Insert elements in map
offcol = 0;
offrow = 0;
for i=1:length(newp)
    aux_map = map; % save the map at its current state (for potential 
    % relocation purposes)
    indel = indp{i};

    % Update index position
    indel(1) = indel(1) + offrow;
    indel(2) = indel(2) + offcol;

    % If the index is in the boundaries (or even further) of the map, then
    % we will have to relocate the elements in map and create a bigger grid
    offrow0 = offrow;
    offcol0 = offcol;
    if indel(1) >= size(map, 1) % last row or more
        nrows = 1 + indel(1) - size(map, 1); % number of additional rows 
        % in the grid (last rows)

        % Map relocation
        map = cell(size(map, 1) + nrows, size(map, 2));
        map(1:size(aux_map, 1), :) = aux_map;
        map((size(aux_map, 1) + 1):end, :) = {[NaN NaN]};
    elseif indel(1) <= 1 % first row or less
        nrows = 2 - indel(1); % number of additional rows in the grid
        % (first rows)
        offrow = offrow + nrows; % rows offset

        % Map relocation
        map = cell(size(map, 1) + nrows, size(map, 2));
        map(1:nrows, :) = {[NaN NaN]};
        map((nrows + 1):end, :) = aux_map;
    end
    
    aux_map = map;
    if indel(2) >= size(map, 2) % last column or more
        ncols = 1 + indel(2) - size(map, 2); % number of additional columns
        % in the grid (last columns)

        % Map relocation
        map = cell(size(map, 1), size(map, 2) + ncols);
        map(:, 1:size(aux_map, 2)) = aux_map;
        map(:, (size(aux_map, 2) + 1):end) = {[NaN NaN]};
    elseif indel(2) <= 1 % first column or less
        ncols = 2 - indel(2); % number of additional columns
        % in the grid (first columns)
        offcol = offcol + ncols; % columns offset

        % Map relocation
        map = cell(size(map, 1), size(map, 2) + ncols);
        map(:, 1:ncols) = {[NaN NaN]};
        map(:, (ncols + 1):end) = aux_map;
    end
    
    % Update the current element index
    indel(1) = indel(1) + (offrow - offrow0);
    indel(2) = indel(2) + (offcol - offcol0);

    % Include elements in map
    map{indel(1), indel(2)} = newp{i};
end
end

% [Future work]: % % Insert elements in tour?
% switch heuristics
%     case 'nearest neighbour'
% 
%         % Find the nearest neighbour in the tour and insert the new
%         % observation before or after this element
%         for i=1:length(newp)
% 
%             % Find the nearest element in the tour
%             mindist = inf;
%             for j=1:length(tour)
%                 dist = norm(tour{j} - newp{i});
%                 if dist < mindist
%                     mindist = dist;
%                     idx     = j;
%                 end
%             end
% 
%             % Evaluate before and after elements distance to decide 
%             % whether the new point should be put before or after its 
%             % nearest neighbour
%             aux_tour = tour;
%             if idx ~= 1 && idx < length(tour)
%                 dist1 = norm(tour{idx - 1} - newp{i});
%                 dist2 = norm(tour{idx + 1} - newp{i});
%                 if dist1 < dist2
%                     tour(idx)     = newp(i);
%                     tour(idx+1:end+1) = aux_tour(idx:length(tour));
%                 else
%                     tour(idx+1)     = newp(i);
%                     tour(idx+2:length(tour)+1) = aux_tour(idx+1: ...
%                         length(tour));
%                 end
% 
%             % Particular cases: when the nearest neighbour is the first
%             % or last element in the tour
%             elseif idx == 1
%                 dist = norm(tour{2} - newp{i});
%                 if dist < mindist
%                     tour(2) = newp(i);
%                     tour(3:length(tour)+1) = aux_tour(2:length(tour));
%                 else
%                     tour(1) = newp(i);
%                     tour(2:length(tour)+1) = aux_tour(:);
%                 end
%             else
%                 dist = norm(tour{idx - 1} - newp{i});
%                 if dist < mindist
%                     tour(idx)     = newp(i);
%                     tour(idx + 1) = aux_tour(idx);
%                 else
%                     tour(idx + 1) = newp(i);
%                 end
%             end
%         end
%     case 'manhattan'
% 
%         % Analyze the tour length when inserting the point at each 
%         % position in the tour
%         for i=1:length(newp)
%             mindelta = inf;
%             for j=1:length(tour)
%                 if j ~= length(tour)
%                     %oldd = norm(tour{j+1} - tour{j});
%                     oldd = manhattan(tour{j+1}, tour{j});
%                     %newd = norm(tour{j+1} - newp{i}) + ...
%                     %    norm(tour{j} - newp{i});
%                     newd = manhattan(tour{j+1}, newp{i}) + ...
%                         manhattan(tour{j}, newp{i});
%                     delta = newd - oldd; % increase of tour length when 
%                     % adding the new insertion point
%                 else
%                     %delta = norm(tour{j} - newp{i});
%                     delta = manhattan(tour{j}, newp{i});
%                 end
%                 if delta < mindelta
%                     mindelta = delta;
%                     idx = j;
%                 end
%             end
% 
%             % Insert element in the position that minimizes the tour
%             % deviation length
%             if idx ~= length(tour)
%                 aux_tour = tour;
%                 tour(idx+1) = newp(i);
%                 tour(idx+2:end+1) = aux_tour(idx+1:end);
%             else
%                 tour(end+1) = newp(i);
%             end
% 
%         end
%     case 'simulated annealing'
% 
% end