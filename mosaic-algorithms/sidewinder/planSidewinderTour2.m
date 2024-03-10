function [grid, origin, itour, grid_topo, tour, dirx, diry, dir1, dir2] = ...
planSidewinderTour2(target, roi, sc, inst, inittime, ovlapx, ovlapy, angle)

% Pre-allocate variables
origin = [0, 0];
% Point camera at ROI's centroid
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));

% Intersect ROI with focal plane
targetArea = topo2inst(roi, cx, cy, target, sc, inst, inittime);

% Closest polygon side to the spacecraft's ground track position (this
% will determine the coverage path)
gt1 = groundtrack(sc, inittime, target);
gt2 = groundtrack(sc, inittime + 500, target);
gt1 = topo2inst(gt1, cx, cy, target, sc, inst, inittime);
gt2 = topo2inst(gt2, cx, cy, target, sc, inst, inittime + 500);
[dir1, dir2] = closestSide2(gt1, gt2, targetArea, angle);

% Build reference tile
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
[~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
fpref.width = width;
fpref.height = height;
fpref.angle = angle;

% Focal plane grid discretization
[grid, dirx, diry] = grid2D(fpref, ovlapx, ovlapy, origin, targetArea);

% Boustrophedon decomposition
itour = boustrophedon(grid, dir1, dir2);

% Grid and tour in topographical coordinates
grid_topo = inst2topo(grid, cx, cy, target, sc, inst, inittime);
tour = inst2topo(itour, cx, cy, target, sc, inst, inittime);
indel = [];
for i=1:numel(tour)
    if isempty(tour{i}), indel = [indel i]; end
end
tour(indel) = [];
end