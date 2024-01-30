function [grid, origin, itour, grid_topo, tour, dirx, diry] = planSidewinderTour2(target, roi, sc, inst, inittime, ovlapx, ovlapy, dir1, dir2, angle)

% Pre-allocate variables
tour = {};
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
grid = [];
origin = [0, 0];
% Point camera at ROI's centroid
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));

% Get camera's FOV boundaries and boresight (when pointing at centroid)
[fovbounds, boresight, rotmat] = instpointing(inst, target, sc, inittime, cx, cy);

% Build focal plane
vertex = cspice_spkpos(sc, inittime, targetframe, 'NONE', target);
point = vertex + fovbounds(:, 1);
plane = cspice_nvp2pl(boresight, point);

% Intersect ROI with focal plane
spoint = zeros(size(roi, 1), 3);
for i=1:size(roi, 1)
    dir = -trgobsvec(roi(i, :), inittime, target, sc);
    [found, spoint(i, :)] = cspice_inrypl(vertex, dir, plane);
    if found == 0
        disp("No intersection");
    end
end

% Build reference tile
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
[~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
fpref.width = width;
fpref.height = height;
fpref.angle = angle;

% Build grid 2D in the focal plane
tArea = zeros(length(spoint), 3);
for i=1:length(spoint)
    vpoint = -(vertex - spoint(i, :)');
    tArea(i, :) = rotmat\vpoint;
end
targetArea = tArea(:, 1:2);

% Focal plane grid discretization
[grid, dirx, diry] = grid2D(fpref, ovlapx, ovlapy, origin, targetArea);

% Boustrophedon decomposition
itour = boustrophedon(grid, dir1, dir2);

% Grid and tour in topographical coordinates
grid_topo = inst2topo(grid, cx, cy, target, sc, inst, inittime);
tour = inst2topo(itour, cx, cy, target, sc, inst, inittime);

end