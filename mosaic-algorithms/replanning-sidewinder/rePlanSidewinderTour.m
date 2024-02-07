function [tour, grid_topo] = rePlanSidewinderTour(target, roi, sc, inst, et, ...
    ovlapx, ovlapy, angle, cx, cy, grid_topo, origin_topo, old_origin_topo)

% Pre-allocate variables
tour = {};
persistent fpref;
persistent dir1;
persistent dir2;
if isempty(origin_topo)
    % Point camera at ROI's centroid
    [origin_topo(1), origin_topo(2)] = centroid(polyshape(roi(:, 1), roi(:, 2)));
end

% Build reference tile (it's always going to be the same in the subsequent
% calls)
if isempty(fpref)
    [~, ~, ~, bounds] = ...
        cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
    [~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
    fpref.width = width;
    fpref.height = height;
    fpref.angle = angle;
end

% Intersect ROI with focal plane
targetArea = topo2inst(roi, cx, cy, target, sc, inst, et);
% Project origin to focal plane
origin = topo2inst(origin_topo, cx, cy, target, sc, inst, et);

% Get sweeping directions (if necessary)
if isempty(dir1)
    % Closest polygon side to the spacecraft's ground track position (this
    % will determine the coverage path)
    gt1 = groundtrack(sc, et, target);
    gt2 = groundtrack(sc, et + 500, target);
    gt1 = topo2inst(gt1, cx, cy, target, sc, inst, et);
    gt2 = topo2inst(gt2, cx, cy, target, sc, inst, et + 500);
    [dir1, dir2] = closestSide2(gt1, gt2, targetArea, angle);
end

if isempty(grid_topo) % First iteration
    % Focal plane grid discretization
    grid = grid2D(fpref, ovlapx, ovlapy, origin, targetArea);
else
    % Replanning Sidewinder: unlike Sidewinder, which only calls
    % planSidewinderTour once, replanning Sidewinder "replans" after
    % every iteration (i.e., footprint). Therefore, in order not to
    % lose the original coverage path direction, a boolean variable is
    % defined to identify when the tour is hopping to a new row/column.
    % If this variable was not put, planSidewinderTour would always
    % start at the right/left or top/bottom side of the grid (depending
    % on the spacecraft's moving direction). So, for example, if we
    % have a horizontal sweep, instead of doing
    % this:
    % ---->-----|
    %           |
    % |---<-----|
    % |
    % |--->------
    % It would be doing this (discontinuities between lines)
    % ----------->
    % >-----------
    % ------------>
    %
    % Likewise, if we have a vertical sweep, instead of doing this:
    % |  ----
    % |  |  |
    % ⌄  ^  ⌄
    % |  |  |
    % |  |  |
    % ----  |
    % It would be doing this (discontinuities between lines)
    % | | |
    % | | |
    % ⌄ ⌄ ⌄
    % | | |
    % | | |
    % | | |
    %
    % The following algorithm identifies if the next element in the
    % tour will imply a change of row/column in the grid
    [i, j] = find(~cellfun('isempty', grid_topo));
    old = []; curr = [];
    for k=1:numel(i)
        if isequal(old_origin_topo, grid_topo{i(k), j(k)})
            old   = [i(k), j(k)];
        elseif isequal(origin_topo, grid_topo{i(k), j(k)})
            curr  = [i(k), j(k)];
        end
        if ~isempty(old) && ~isempty(curr)
            break;
        end
    end
    
    if isequal(dir1, 'east') || isequal(dir1, 'west') % vertical sweep
        if curr(2) ~= old(2)
            % direction must be the contrary to the current one
            if isequal(dir2, 'north'), dir2 = 'south';
            else, dir2 = 'north'; end
        end
    elseif isequal(dir1, 'north') || isequal(dir1, 'south')
        if curr(1) ~= old(1) % change of row between replans
            % next row's coverage path's
            if isequal(dir2, 'west'), dir2 = 'east';
            else, dir2 = 'west'; end
        end
    end

    % Replanning sidewinder: optimize grid origin in order to avoid
    % potential taboo tiles
    grid = optimizeGridOrigin2(origin, fpref, ovlapx, ovlapy, ...
          targetArea, dir1, dir2);
    if isempty(grid), return; end

    % flag = true;
    % gamma0 = origin';
    % grid0  = grid;
    % while flag
    %     grid = optimizeGridOrigin2(gamma0, fpref, ovlapx, ovlapy, ...
    %         targetArea, dir1, dir2);
    %     if ~isempty(grid) 
    %         flag = false;
    %     else 
    %         gamma = gamma0;
    %         count = 0;
    %         delta = [0, 0];
    %         deltax = 0.1*fpref.width; % displacement value in x direction
    %         deltay = 0.1*fpref.height; % displacement value in y direction
    %         while isempty(grid) && count ~= 10
    %             count = count + 1;
    %             [cx, cy] = centroid(polyshape(targetArea(:, 1), targetArea(:, 2)));
    %             dirc = [cx - gamma0(1), cy - gamma0(2)];
    %             dirc = dirc/norm(dirc);
    %             delta = delta + [deltax*dirc(1), deltay*dirc(2)];
    %             gamma = gamma + delta;
    %             grid = grid2D(fpref, ovlapx, ovlapy, gamma, targetArea);
    %         end
    % 
    %         if ~isempty(grid)
    %             % Find which position does gamma occupy in this grid
    %             for i=1:numel(grid)
    %                 if ~isempty(grid{i}) && norm(grid{i} - gamma') < 1e-3
    %                     [ind_row, ind_col] = ind2sub(size(grid), i);
    %                     break;
    %                 end
    %             end
    %             grid0{ind_row, ind_col} = grid{ind_row, ind_col};
    %             grid = grid0;
    %             flag = false;
    %         else
    %             return;
    %         end
    %     end
        
end

% Boustrophedon decomposition
itour = boustrophedon(grid, dir1, dir2);

% Transform coordinates
grid_topo = inst2topo(grid, cx, cy, target, sc, inst, et);
tour = inst2topo(itour, cx, cy, target, sc, inst, et);
indel = [];
for i=1:numel(tour)
    if isempty(tour{i}), indel = [indel i]; end
end
tour(indel) = [];
end