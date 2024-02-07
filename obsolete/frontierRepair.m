function [A, coverage, fpList] = frontierRepair(startTime, endTime, ...
    tobs, inst, sc, target, roi, olapx, olapy, ax, c, video)
% [Description]
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% Last Rev.     07/2023
% 
% Usage:        [A, coverage, fpList] = frontierRepair(startTime, endTime, ...
%               tobs, inst, sc, target, roi, olapx, olapy, ax, c, video)
%
% Inputs:
%   > startTime:    start time of the planning horizon, in TDB seconds past
%                   J2000 epoch
%   > endTime:      end time of the planning horizon, in TBD seconds past
%                   J2000 epoch
%   > tobs:         observation time, i.e. the minimum time that the 
%                   instrument needs to perform an observation, in seconds
%   > inst:         string name of the instrument
%   > sc:           string name of the spacecraft
%   > target:       string name of the target body
%   > roi:          matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D. 
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > ax:
%   > c:
%   > video:
% 
% Outputs:
%   > A:            cell matrix of the successive instrument observations,
%                   sorted in chronological order.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > coverage:
%   > fpList:
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables
A = {}; % List of observations (successive boresight ground track position)
theta = 0; % temppppppp
[~, targetFrame, ~] = cspice_cnmfrm(target); % body-fixed frame
fpList = [];
if ~isempty(ax)
    mapplot = 1;
else
    mapplot = 0;
end

% Previous anti-meridian intersection check...
ind = find(diff(sort(roi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind)
    [x, y] = amsplit(roi(:, 1), roi(:, 2));
    roi(roi(:, 1) < 0, 1) = roi(roi(:, 1) < 0, 1) + 360;
    [roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2));
else
    x = roi(:,1); y = roi(:,2);
end

% Define target area as a polygon
poly1 = polyshape(x, y);
roiarea = area(poly1); % surface area enclosed by the roi
%roiarea = polysurfarea([x, y], target); 

% Animation of coverage map
if ~isempty(video)
    open(video);
    writeVideo(video,getframe(gcf));
end

%% Frontier Repair algorithm
% The first time iteration is the starting time in the planning horizon
t = startTime;

% Boolean that defines when to stop covering the target area
exit = false;

while ~exit && t < endTime 

    % Initial 2D grid layout discretization: the instrument's FOV is going
    % to be projected onto the uncovered area's centroid and the resulting
    % footprint shape is used to set the grid spatial resolution
    [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1),roi(:,2)));
    fprintc = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
        theta);   % centroid footprint
    % If the footprint contains the limb, the centroid of the footprint
    % might not coincide with the camera's boresight projection
    if fprintc.limb
        gamma(1) = fprintc.clon; gamma(2) = fprintc.clat;
    end
    
    % Initialize struct that saves footprints (sub-structs)
    if t == startTime
        fpList = struct([]);
        for fn = fieldnames(fprintc)'
            fpList(1).(fn{1}) = [];
        end
    end

    if isempty(fprintc.bvertices)
        disp("Region of interest not visible from the instrument")
        exit = true;
        continue % the footprint is empty because the roi (in this case,
        % more specifically, the center point of the roi area) is not
        % visible from the instrument FOV. Therefore, the function is
        % exited. This may not be completely correct because the roi area
        % could be large enough so that its centroid is not visible but any
        % other region inside the area is. This is left for future work.
    end

    % Closest polygon side to the spacecraft's ground track position (this
    % will determine the coverage path in planSidewinderTour)
    cside = closestSide(target, sc, t, roi, fprintc.angle);

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    [tour, ~, sweep] = planSidewinderTour(target, sc, t, roi, fprintc,...
        olapx, olapy, cside, {gamma}, {});
    if isempty(tour)
        exit = true;
        continue
    end
    clear planSidewinderTour;
    [grid, dirx, diry] = grid2D(fprintc, olapx, olapy, gamma, roi);

    % Grid map: this map encloses the discretized grid. First and last rows
    % and columns are NaN, as well as the empty values in the grid (points
    % outside of the ROI or excluded from observation bc their associated
    % allocated cell's coverage is too small)
    map = cell(size(grid,1) + 2, size(grid,2) + 2);
    for i=1:size(map,1)
        for j=1:size(map,2)
            if i == 1 || j == 1 || i == size(map,1) || ...
                    j == size(map,2) || any(isempty(grid{i-1, j-1}))
                map{i,j} = [NaN NaN];
            else
                map{i, j} = grid{i-1, j-1}';
            end
        end
    end

    a = []; % initialize
    while ~isempty(tour)

        gamma  = tour{1}; % current observation point
        gamma0 = a; % previous observation point

        % Compute the observation point's footprint
        currfp = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
            theta);
        % If the footprint contains the limb, the centroid of the footprint
        % might not coincide with the camera's boresight projection
        if currfp.limb
            gamma(1) = currfp.clon; gamma(2) = currfp.clat;
        end

        % If the footprint is no longer visible, then move on to the next
        % planned observation
        if isempty(currfp.bvertices)
            tour(1) = [];
            continue;
        end

        % Update previous grid with the new tile reference (footprint),
        % looking for new potential tiles and/or disposable ones
        [N, Nind, X, map, tour, sweep] = updateGrid(roi, tour, gamma, ...
            gamma0, cside, sweep, map, olapx, olapy, dirx, diry, currfp);

        % Remove disposable tiles
        [tour, map] = removeTiles(tour, map, X);

        % Insert new tiles
        [tour, map] = insertTiles(tour, map, N, Nind);

        if isempty(tour)
            continue;
        end

        % Update tour
        a = tour{1}; % observation
        tour(1) = []; % delete this observation from the planned tour

        % Check a.m. intercept...
        if a(1) > 180, a(1) = a(1) - 360; end

        % If the current observation point has been deleted, then we have
        % to re-calculate the current footprint
        if norm(a - gamma) > 1e-5
            currfp = footprint(a(1), a(2), t, inst, sc, target, ...
            theta);
        end
        
        % Subtract the corresponding area from the target polygon
        A{end + 1} = a; % add it in the list of planned observations
        poly2 = polyshape(currfp.bvertices); % create footprint polygon
        poly1 = subtract(poly1, poly2); % update uncovered area
        roi = poly1.Vertices;

        % Save footprint struct
        fpList(end + 1) = currfp;

        %% Plots
        if mapplot
            % Footprint plot in Figure 2
            plot(ax, poly2, 'FaceColor', 'b', 'EdgeColor', ...
                'b', 'linewidth', 1, 'FaceAlpha', 0.2)
            if length(A) > 1
                if abs(A{end-1}(1) - A{end}(1)) <= 180 % no coverage ...
                    % path - a.m. intercept
                    plot(ax, [A{end-1}(1) A{end}(1)], [A{end-1}(2) A{end}(2)],...
                        'w-', 'linewidth', 1)
                end
            end
            drawnow
        end

        % Figure showing the footprint projection onto the body surface
        %footprint3Dprojection(currfp, videosave)

        % Spacecraft ground track position in Figure 2
        if mapplot
            sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
                t, targetFrame, 'NONE', sc);
            [~, sclon, sclat] = cspice_reclat(sctrack);
            scatter(ax, sclon*cspice_dpr, sclat*cspice_dpr, 8, c,...
                "filled")
        end

        % Save Figure 2 frame in the animation
        if ~isempty(video), writeVideo(v2, getframe(gcf)); end

        % New time iteration
        t = t + tobs;
        %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work
    end

    %% Future work: automated scheduling (in-situ)
    % % Stop criteria: if the surface of the remaining uncovered roi area is
    % % smaller than half of the last footprint size, then it is not worth it
    % % to start the tour (for the uncovered roi area) again
    % if area(polyshape(poly1.Vertices(:, 1), poly1.Vertices(:, 2))) < ...
    %         0.2*area(polyshape(lastfp.bvertices(:, 1), ...
    %         lastfp.bvertices(:, 2)))
    %     exit = true;
    % else
    %     roi = poly1.Vertices;
    % end

    % For now, the stop criteria is the end of the tour, re-starts are not
    % optimal for the purposes of the scheduling problem
    exit = true;
end

% Remove first element of fplist (it was just to set the struct fields)
if ~isempty(fpList)
    fpList(1) = [];
end

% ROI coverage percentage
%coverage = (roiarea - polysurfarea(poly1.Vertices, target))/roiarea;
coverage = (roiarea - area(poly1)) / roiarea;

% End animations
%if videosave, close(v2); end
end