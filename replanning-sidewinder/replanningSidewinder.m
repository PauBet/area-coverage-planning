function [A, coverage, fpList] = replanningSidewinder(startTime, endTime, ...
    tobs, inst, sc, target, roi, olapx, olapy, ax, c, videosave)
% Optimization of the Sidewinder algorithm, adapted from [1]. 
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        A = replanningSidewinder(startTime, endTime, tobs,...
%                   inst, sc, target, roi, olapx, olapy, videosave)
%
% Inputs:
%   > startTime:    start time of the planning horizon, in TDB seconds past
%                   J2000 epoch
%   > endTime:      end time of the planning horizon, in TBD seconds past
%                   J200 epoch
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
%   > videosave:
% 
% Outputs:
%   > A:            cell matrix of the successive instrument observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

%%
% Pre-allocate variables
A = {}; % List of observations (successive boresight ground track position)
grid = {}; % Sorted ROI's grid discretization (Boustrophedon decomposition)
exit = false; % Boolean that defines when to stop covering the target area
theta = 0; % temppppppp
[~, targetFrame, ~] = cspice_cnmfrm(target); % body-fixed frame
start = true; % boolean that indicates the start of the algorithm

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
poly1 = polyshape(x, y); % remaining area (may not be the same as unroi
% because sometimes the footprints may be blank and, yet, the "allocated"
% area for the footprint must be removed to avoid gridlock)
roiarea = area(poly1); % surface area enclosed by the roi
unroi = poly1; % total uncovered area

% Animation of coverage map
if videosave
    v2 = VideoWriter('topography_map_rsidewinder');
    v2.FrameRate = 2;
    open(v2);
    writeVideo(v2,getframe(gcf));
end

%% Replanning Sidewinder algorithm

% The first time iteration is the starting time in the planning horizon
t = startTime;
% Initial 2D grid layout discretization: the instrument's FOV is
% going to be projected onto the uncovered area's centroid and the
% resulting footprint shape is used to set the grid spatial
% resolution
[gamma(1), gamma(2)] = centroid(polyshape(roi(:,1), roi(:,2)));
tour{1} = gamma;
lastfp  = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
    theta);  % reference footprint at gamma

% Initialize struct that saves footprints (sub-structs)
fpList = struct([]);
for fn = fieldnames(lastfp)'
   fpList(1).(fn{1}) = [];
end

% Closest polygon side to the spacecraft's ground track position
cside = closestSide(target, sc, t, roi);

while ~iszero(roi) && t <= endTime && ~exit

    %fprint0 = footprint(gamma(1), gamma(2), t, inst, sc, ...
    %    target, theta);  % reference footprint at gamma
    fprint0 = lastfp;

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    if length(tour) > 2 || start
        [tour, grid] = planSidewinderTour(target, sc, t, roi, fprint0,...
            olapx, olapy, cside, tour, grid);
        start = false;
    else
        tour(1) = []; % if it is the last element of the tour, there is no
        % need to replan. In case there should have been a necessity 
        % because the footprint still leaves a significant portion of roi
        % uncovered (worst case scenario), then it will re-start with the 
        % uncovered area. This has been done like this because, sometimes,
        % roi may be divided in more than one polygon, and so the gamma
        % gets isolated and the flood fill algorithm cannot discretize the
        % entire region of uncovered roi (it simply does not find gamma's
        % neighbors because they're too separated from it, think about the
        % last portion of a star tip, for example). In order to avoid this
        % problem, planSidewinderTour has identified these cases as
        % length(tour) == 1, which could interfere with the case where
        % length(tour) == 1 because it's simply the last footprint of the
        % discretized area
    end

    if ~isempty(tour)
        % Compute the footprint of each point in the tour successively and
        % subtract the corresponding area from the target polygon
        a = tour{1}; % observation

        % Check a.m. intercept...
        if a(1) > 180, a(1) = a(1) - 360; end

        % Compute the observation's footprint
        fprinti = footprint(a(1), a(2), t, inst, sc, target, ...
            theta);

        if ~isempty(fprinti.bvertices)
            A{end + 1} = a; % add it in the list of planned observations
            poly2 = polyshape(fprinti.bvertices); % create footprint 
            % polygon
            poly1 = subtract(poly1, poly2); % update remaining area
            unroi = subtract(unroi, poly2); % update total uncovered area

            % Save footprint struct
            fpList(end + 1) = fprinti;

            %% Plots
            % Footprint plot in Figure 2
            plot(ax, poly2, 'FaceColor', 'b', 'EdgeColor', ...
                'b', 'linewidth', 1, 'FaceAlpha', 0.2)
            if length(A) > 1
                if abs(A{end-1}(1) - A{end}(1)) <= 180 % no coverage path -
                    % a.m. intercept
                    plot(ax, [A{end-1}(1) A{end}(1)], [A{end-1}(2) ...
                        A{end}(2)], 'w-', 'linewidth', 1)
                end
            end
            drawnow

            % Figure showing the footprint projection onto the body surface
            %footprint3Dprojection(fprinti, videosave)

            % Spacecraft ground track position in Figure 2
            groundTrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
                t, targetFrame, 'NONE', sc);
            [~, sclon, sclat] = cspice_reclat(groundTrack);
            plot(ax, sclon*cspice_dpr, sclat*cspice_dpr, 'y.', ...
                'MarkerSize', 6)

            % Save Figure 2 frame in the animation
            if videosave, writeVideo(v2, getframe(gcf)); end

            % New time iteration
            t = t + tobs;
            %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work

            lastfp = fprinti;
            roi    = poly1.Vertices;
        else
            disp("Footprint not visible from the instrument")

            % Erase the allocated area of the roi (to prevent the algorithm
            % to go backwards)
            faux = lastfp;
            faux.bvertices(:, 1) = faux.bvertices(:, 1) - faux.olon + a(1);
            faux.bvertices(:, 2) = faux.bvertices(:, 2) - faux.olat + a(2);
            faux.olon = a(1); faux.olat = a(2);

            poly2 = polyshape(faux.bvertices); % create footprint
            % polygon
            poly1 = subtract(poly1, poly2); % update uncovered area
            roi    = poly1.Vertices;
        end
    else
        % Stop criteria: if the surface of the remaining uncovered roi area
        % is smaller than half of the last footprint size, then it is not
        % worth it to start the tour (for the uncovered roi area) again
        if area(polyshape(poly1.Vertices(:, 1), poly1.Vertices(:, 2))) < ...
                0.2*area(polyshape(lastfp.bvertices(:, 1), ...
                lastfp.bvertices(:, 2)))
            exit = true;
        else
            roi   = poly1.Vertices;
            [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1), roi(:,2)));
            lastfp  = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
                theta);  % reference footprint at gamma
            tour{1} = gamma;
            start = true;
            
            % COMMENT THIS
            polyfp = polyshape(lastfp.bvertices(:, 1), ...
                lastfp.bvertices(:, 2));
            polyout = intersect(polyfp, poly1);
            if area(polyout) < 0.2*area(polyfp)
                exit = true;
            end
        end
    end
end

% Remove first element of fplist (it was just to set the struct fields)
fpList(1) = [];

% ROI coverage percentage
%coverage = (roiarea - polysurfarea(poly1.Vertices, target))/roiarea;
coverage = (roiarea - area(unroi)) / roiarea;

% End animations
if videosave, close(v2); end
end