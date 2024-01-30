function [tour, grid_topo] = rePlanSidewinderTour(target, roi, sc, inst, et, ...
    ovlapx, ovlapy, angle, indir1, indir2, grid_topo, origin_topo, old_origin_topo)

% Pre-allocate variables
tour = {};
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
persistent fpref;
persistent dir1;
persistent dir2;
if isempty(dir1)
    dir1 = indir1; dir2 = indir2;
end
if isempty(origin_topo)
    % Point camera at ROI's centroid
    [origin_topo(1), origin_topo(2)] = centroid(polyshape(roi(:, 1), roi(:, 2)));
end
method = 'ELLIPSOID'; % assumption: ray intercept function is going to
% model the target body as a tri-axial ellipsoid

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

% Get camera's FOV boundaries and boresight (when pointing at origin)
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
[fovbounds, boresight, rotmat] = instpointing(inst, target, sc, et, cx, cy);

% Build focal plane
vertex = cspice_spkpos(sc, et, targetframe, 'NONE', target);
point = vertex + fovbounds(:, 1);
plane = cspice_nvp2pl(boresight, point);

% Intersect ROI with focal plane
spoint = zeros(size(roi, 1), 3);
for i=1:size(roi, 1)
    if ~isnan(roi(i, :))
        dir = -trgobsvec(roi(i, :), et, target, sc);
        [found, spoint(i, :)] = cspice_inrypl(vertex, dir, plane);
        if found == 0
            disp("No intersection");
        end
    else
        spoint(i, :) = nan(1, 3);
    end
end

% Build grid 2D in the focal plane
tArea = zeros(length(spoint), 3);
for i=1:length(spoint)
    if ~isnan(spoint(i, :))
        vpoint = -(vertex - spoint(i, :)');
        tArea(i, :) = rotmat\vpoint;
    else
        tArea(i, :) = nan(1, 3);
    end
end
targetArea = tArea(:, 1:2);

% Intersect origin with focal plane
dir = -trgobsvec(origin_topo, et, target, sc);
[found, aux] = cspice_inrypl(vertex, dir, plane);
if found == 0
    % Pending work...
    disp("Origin not visible");
    return;
end
vpoint = -(vertex - aux);
aux = rotmat\vpoint;
origin = aux(1:2, :);

if isempty(grid_topo) % First iteration
    % Focal plane grid discretization
    grid = grid2D(fpref, ovlapx, ovlapy, origin', targetArea);
else
    % Intersect old_origin with focal plane
    dir = -trgobsvec(old_origin_topo, et, target, sc);
    [found, aux] = cspice_inrypl(vertex, dir, plane);
    if found == 0
        % This is a special case... pending work
        disp("Old origin not visible");
        return;
    end
    vpoint = -(vertex - aux);
    aux = rotmat\vpoint;
    old_origin = aux(1:2, :);

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
    grid = optimizeGridOrigin2(origin', fpref, ovlapx, ovlapy, ...
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
% 
% figure
% plot(polyshape(targetArea(:, 1), targetArea(:, 2)))
% hold on; axis equal;
% plot(targetArea(:, 1), targetArea(:, 2), 'b*')
% for i=1:size(grid, 1)
%     for j=1:size(grid, 2)
%         sp = grid{i, j};
%         if ~isempty(sp)
%             plot(sp(1), sp(2), '^')
%         end
%     end
% end

grid_topo = inst2topo(grid, rotmat, target, et, sc);

% Boustrophedon decomposition
itour = boustrophedon(grid, dir1, dir2);

% Intersect each tile with the target's surface
count = 0;
for i=1:length(itour)
    sp = itour{i};
    if ~isempty(sp)
        p = zeros(3, 1);
        p(1:2) = sp; p(3) = 1;
        p_body = rotmat*p;
        [xpoint, ~, ~, found] = cspice_sincpt(method, target, et,...
            targetframe, 'NONE', sc, targetframe, p_body);
        if found
            count = count + 1;
            [~, lon, lat] = cspice_reclat(xpoint);
            tour{count} = [lon*cspice_dpr, lat*cspice_dpr];
        end
    end
end

end