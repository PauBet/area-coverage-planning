function [A, coverage, fpList] = replanningSidewinder(startTime, endTime, ...
    tobs, inst, sc, target, roi, olapx, olapy, ax, video)
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
%   > olapx:        grid footprint overlap in the x direction,
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction,
%                   in percentage (height)
%   > ax:           axes of the figure (map) where the subsequent 
%                   observations are going to be projected
%   > video:        video object where the evolution of the figure (ax)
%                   over the observation plan is going to be animated
% 
% Outputs:
%   > A:            cell matrix of the successive instrument observations.
%                   Each observation is defined by the instrument 
%                   boresight projection onto the body surface, in 
%                   latitudinal coordinates [lon lat], in deg
%   > coverage:     total observation coverage (coveredarea/totalarea)
%   > fpList:       struct that containts the list of footprints 
%                   (observations)
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate and define variables
A = {}; % List of observations (successive boresight ground track position)
grid = {}; % Sorted ROI's grid discretization (Boustrophedon decomposition)
exit = false; % Boolean that defines when to stop covering the target area
theta = 0; % temppppppp
[~, targetFrame, ~] = cspice_cnmfrm(target); % body-fixed frame
fpList = struct([]);
coverage = 0;

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
if ~isempty(video)
    open(video);
    writeVideo(video,getframe(gcf));
end

%% Replanning Sidewinder algorithm

% The first time iteration is the starting time in the planning horizon
t = startTime;
% Initial 2D grid layout discretization: the instrument's FOV is
% going to be projected onto the uncovered area's centroid and the
% resulting footprint shape is used to set the grid spatial
% resolution
[gamma(1), gamma(2)] = centroid(polyshape(roi(:,1), roi(:,2)));
tour{1} = gamma; % initial tour seed
fprint0  = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
    theta);
% If the footprint contains the limb, the centroid of the footprint
% might not coincide with the camera's boresight projection
% if fprint0.limb
%     gamma(1) = fprint0.clon; gamma(2) = fprint0.clat;
% end

if isempty(fprint0.bvertices)
    disp("ROI not visible... Exiting algorithm")
    return;
end

% Initialize struct that saves footprints (sub-structs)
for fn = fieldnames(fprint0)'
   fpList(1).(fn{1}) = [];
end

% Closest polygon side to the spacecraft's ground track position (this 
% sets the starting point of the tour)
cside = closestSide(target, sc, t, roi, fprint0.angle);

% Start using the footprint that corresponds to the first point in the
% tour
[tour, grid] = planSidewinderTour(target, sc, t, roi, fprint0,...
    olapx, olapy, cside, tour, grid);
currfp  = footprint(tour{1}(1), tour{1}(2), t, inst, sc, target, ...
    theta);  % reference footprint at gamma
if currfp.limb
    gamma(1) = currfp.clon; gamma(2) = currfp.clat;
end
clear planSidewinderTour;
tour = [];

while ~iszero(roi) && t <= endTime && ~exit

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    %if length(tour) > 2 || start
    if length(tour) > 1
        [tour, grid] = planSidewinderTour(target, sc, t, roi, currfp,...
            olapx, olapy, cside, tour, grid);
    else
        [tour, grid] = planSidewinderTour(target, sc, t, roi, currfp,...
            olapx, olapy, cside, {gamma}, []);
        % tour(1) = []; % last element of the tour was already observed
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
            if ~isempty(ax)
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
            end

            % Figure showing the footprint projection onto the body surface
            % (auxiliary figure)
            %footprint3Dprojection(fprinti, 0)

            % Spacecraft ground track position in Figure 2
            if ~isempty(ax)
                groundTrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
                    t, targetFrame, 'NONE', sc);
                [~, sclon, sclat] = cspice_reclat(groundTrack);
                plot(ax, sclon*cspice_dpr, sclat*cspice_dpr, 'y.', ...
                    'MarkerSize', 6)
            end

            % Save Figure 2 frame in the animation
            if ~isempty(video), writeVideo(video, getframe(gcf)); end

            % New time iteration
            t = t + tobs;
            %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work

            currfp = fprinti;
            if length(tour) > 1
                nextfp = footprint(tour{2}(1), tour{2}(2), t, inst, sc, target, ...
                    theta);
                if ~isempty(nextfp.bvertices)
                    currfp = nextfp;
                    if currfp.limb
                        % When the footprint contains the limb, it may
                        % happen that the boresight projection of the
                        % camera does not coincide with the centroid of the
                        % footprint, which is the reference point to
                        % create the grid discretization. For this reason,
                        % we shall change the boresight projection for the
                        % centroid, in order to get a better discretization
                        % of the ROI
                        gamma_old(1) = tour{2}(1); gamma_old(2) = tour{2}(2);
                        tour{2}(1) = currfp.clon; tour{2}(2) = currfp.clat;
                        [i, j] = find(~cellfun('isempty', grid));
                        gidx = [];
                        for k=1:numel(i)
                            if isequal(gamma_old', grid{i(k), j(k)})
                                gidx   = [i(k), j(k)];
                            end
                            if ~isempty(gidx)
                                break;
                            end
                        end
                        grid{gidx(1), gidx(2)} = tour{2}';
                    end
                end
            end
            roi    = poly1.Vertices;
        else
            disp("Footprint not visible from the instrument")

            % Erase the allocated area of the roi (to prevent the algorithm
            % going backwards)
            faux = currfp;
            faux.bvertices(:, 1) = faux.bvertices(:, 1) - faux.olon + a(1);
            faux.bvertices(:, 2) = faux.bvertices(:, 2) - faux.olat + a(2);
            faux.olon = a(1); faux.olat = a(2);

            poly2 = polyshape(faux.bvertices); % create footprint
            % polygon
            poly1 = subtract(poly1, poly2); % update uncovered area
            roi   = poly1.Vertices;
        end
    else
        % Stop criteria: if the surface of the remaining uncovered roi area
        % is smaller than half of the last footprint size, then it is not
        % worth it to start the tour (for the uncovered roi area) again
        if area(polyshape(poly1.Vertices(:, 1), poly1.Vertices(:, 2))) < ...
                0.2*area(polyshape(currfp.bvertices(:, 1), ...
                currfp.bvertices(:, 2)))
            exit = true;
        else
            % COMMENT THIS
            polyfp = polyshape(currfp.bvertices(:, 1), ...
                currfp.bvertices(:, 2));
            polyout = intersect(polyfp, poly1);
            if area(polyout) < 0.2*area(polyfp)
                exit = true;
                continue
            end

            roi   = poly1.Vertices;
            [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1), roi(:,2)));
            currfp  = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
                theta);  % reference footprint at gamma
            tour{1} = gamma;
        end
    end
end

% Remove first element of fplist (it was just to set the struct fields)
fpList(1) = [];

% ROI coverage percentage
coverage = (roiarea - area(unroi)) / roiarea;

% End animations
%if ~isempty(video), close(v2); end
end