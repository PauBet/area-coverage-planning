function [tour, roi] = get_tour_first(inittime, inst, sc, ...
    target, inroi)

    ovlapx = 20;
    ovlapy = 20;
    
    % Check ROI visible area from spacecraft
    vsbroi = visibleroi(inroi, inittime, target, sc);
    roi = interppolygon(vsbroi); % interpolate polygon vertices (for improved 
    % accuracy)

    % Previous anti-meridian intersection check...
    ind = find(diff(sort(inroi(:, 1))) >= 180, 1); % find the discontinuity index
    if ~isempty(ind)
        roi = inroi;
        roi(roi(:, 1) < 0, 1) = roi(roi(:, 1) < 0, 1) + 360;
        [roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2));
    end

    % Initial 2D grid layout discretization: the instrument's FOV is going
    % to be projected onto the uncovered area's centroid and the resulting
    % footprint shape is used to set the grid spatial resolution

    [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1),roi(:,2)));
    fprintc = footprint(inittime, inst, sc, target, 'lowres', ...
        gamma(1), gamma(2), 1);
    
    angle = fprintc.angle;

    % Pre-allocate variables
    A = {}; % List of observations (successive boresight ground track position)
    fpList = [];

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
    [dir1, dir2] = closestSide(gt1, gt2, targetArea, -angle);
    
    % Build reference tile
    [~, ~, ~, bounds] = ...
        cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
    [~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
    fpref.width = width;
    fpref.height = height;
    fpref.angle = -angle;

    % Focal plane grid discretization
    [grid, ~, ~] = grid2D(fpref, ovlapx, ovlapy, origin, targetArea);

    % Boustrophedon decomposition
    itour = boustrophedon(grid, dir1, dir2);
    
    % Tour in topographical coordinates
    tour = inst2topo(itour, cx, cy, target, sc, inst, inittime);

end
