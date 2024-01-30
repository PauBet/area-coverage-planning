function [A, coverage, fpList, makespan, nfp] = sidewinder(startTime, endTime, tobs,...
    inst, sc, target, roi, olapx, olapy, ax, video, express)
% Target-fixed Boustrophedon decomposition, adapted from [1].
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        A = sidewinder(startTime, endTime, tobs, instName, ...
%               scName, targetName, vertices, method)
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
%                   vertex points are expressed in 2D, in latitudinal 
%                   coordinates [º]
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > videosave:
% 
% Outputs:
%   > A:            cell matrix of the successive instrument observations,
%                   sorted in chronological order.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

%
% Pre-allocate variables
A = {}; % List of observations (successive boresight ground track position)
theta = 0; % temppppppp
[~, targetFrame, ~] = cspice_cnmfrm(target); % body-fixed frame
fpList = [];
nfp = 0;
makespan = 0;
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

%% Sidewinder algorithm

% The first time iteration is the starting time in the planning horizon
t = startTime;

% Boolean that defines when to stop covering the target area
exit = false;

while ~exit && t < endTime 

    % Initial 2D grid layout discretization: the instrument's FOV is going
    % to be projected onto the uncovered area's centroid and the resulting
    % footprint shape is used to set the grid spatial resolution
    [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1),roi(:,2)));
    if ~express
        fprintc = footprint(t, inst, sc, target, 'highres', ...
            gamma(1), gamma(2)); % centroid footprint
    else
        fprintc = footprint(t, inst, sc, target, 'lowres', ...
            gamma(1), gamma(2)); % centroid footprint
    end
    % If the footprint contains the limb, the centroid of the footprint
    % might not coincide with the camera's boresight projection
    if ~isequal(fprintc.limb, 'none')
        gamma(1) = fprintc.clon; gamma(2) = fprintc.clat;
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

    % Initialize struct that saves footprints (sub-structs)
    if t == startTime
        fpList = struct([]);
        for fn = fieldnames(fprintc)'
            fpList(1).(fn{1}) = [];
        end
    end

    % If the footprint is larger than the ROI, stop
    xq = roi(:, 1); yq = roi(:, 2);
    xv = fprintc.bvertices(:, 1); yv = fprintc.bvertices(:, 2);
    if any(inpolygon(xq, yq, xv, yv), 'all')
        fpList(end+1) = fprintc;
        tour = {gamma};
        disp("FOV projection is larger than ROI surface")
        exit = true;
        continue 
    end

    % Closest polygon side to the spacecraft's ground track position (this
    % will determine the coverage path in planSidewinderTour)
    cside = closestSide(target, sc, t, roi, fprintc.angle);

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    tour = planSidewinderTour(target, sc, t, roi, fprintc,...
        olapx, olapy, cside, {gamma}, {});

    if isempty(tour)
        exit = true;
        continue
    end

    if ~express
        while ~isempty(tour)
            % Compute the footprint of each point in the tour successively and
            % subtract the corresponding area from the target polygon
            a = tour{1}; % observation
            tour(1) = []; % delete this observation from the planned tour

            % Check a.m. intercept...
            if a(1) > 180, a(1) = a(1) - 360; end

            % Compute the observation's footprint
            fprinti = footprint(t, inst, sc, target, 'highres', a(1), a(2));

            if ~isempty(fprinti.bvertices)
                A{end + 1} = a; % add it in the list of planned observations
                poly2 = polyshape(fprinti.bvertices); % create footprint polygon
                poly1 = subtract(poly1, poly2); % update uncovered area

                % Save footprint struct
                fpList(end + 1) = fprinti;

                %% Plots
                if mapplot
                    % Footprint plot in Figure 2
                    plot(ax, poly2, 'FaceColor', [0.93, 0.69, 0.13], 'EdgeColor', ...
                        [1, 1, 0.698], 'linewidth', 1, 'FaceAlpha', 0.2)
                    if length(A) > 1
                        if abs(A{end-1}(1) - A{end}(1)) <= 180 % no coverage ...
                            % path - a.m. intercept
                            % plot(ax, [A{end-1}(1) A{end}(1)], [A{end-1}(2) A{end}(2)],...
                            %     'w-', 'linewidth', 1)
                            plot(ax, [A{end-1}(1) A{end}(1)], [A{end-1}(2) A{end}(2)],...
                                'Color', [0.64, 0.08, 0.18], 'linewidth', 1.5)
                        end
                    else
                        scatter(ax, A{1}(1), A{1}(2), 50, 'm', "filled")
                    end
                    drawnow
                end

                % Figure showing the footprint projection onto the body surface
                %footprint3Dprojection(fprinti, 0)

                % Spacecraft ground track position in Figure 2
                if mapplot
                    sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
                        t, targetFrame, 'NONE', sc);
                    [~, sclon, sclat] = cspice_reclat(sctrack);
                    if length(A) > 1
                        scatter(ax, sclon*cspice_dpr, sclat*cspice_dpr, 8, 'c', ...
                            "filled")
                    else
                        scatter(ax, sclon*cspice_dpr, sclat*cspice_dpr, 50, 'm', ...
                            "filled")
                    end
                end

                % Save Figure 2 frame in the animation
                if ~isempty(video), writeVideo(video, getframe(gcf)); end

                % New time iteration
                t = t + tobs;
                %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work

                % Future work: automated scheduling
                % lastfp = fprinti;
            end
        end
    end
    
    %% Future work: automated scheduling (in-situ)
%     % Stop criteria: if the surface of the remaining uncovered roi area is
%     % smaller than half of the last footprint size, then it is not worth it
%     % to start the tour (for the uncovered roi area) again
%     if area(polyshape(poly1.Vertices(:, 1), poly1.Vertices(:, 2))) < ...
%             0.2*area(polyshape(lastfp.bvertices(:, 1), ...
%             lastfp.bvertices(:, 2)))
%         exit = true;
%     else
%         polyfp = polyshape(lastfp.bvertices(:, 1), ...
%             lastfp.bvertices(:, 2));
%         polyout = intersect(polyfp, poly1);
%         if area(polyout) < 0.2*area(polyfp)
%             exit = true;
%             continue
%         end
% 
%         roi   = poly1.Vertices;
%         [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1), roi(:,2)));
%         lastfp  = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
%             theta);  % reference footprint at gamma
%     end
    
    % For now, the stop criteria is the end of the tour, re-starts are not
    % optimal for the purposes of the scheduling problem
    exit = true;
end

if ~express
    % Remove first element of fplist (it was just to set the struct fields)
    if ~isempty(fpList)
        fpList(1) = [];
    else
        return;
    end

    % Make-span
    nfp = length(fpList);
    if length(fpList) > 1
        makespan = fpList(end).t - fpList(1).t;
    else
        makespan = tobs;
        poly1 = polyshape(fpList.bvertices(:, 1), fpList.bvertices(:, 2));
    end

    % ROI coverage percentage
    %coverage = (roiarea - polysurfarea(poly1.Vertices, target))/roiarea;
    coverage = (roiarea - area(poly1)) / roiarea;
else
    A = tour;
    makespan = tobs*length(A);
    nfp = length(A);
    fpList = [];
    coverage = [];
end

% End animations
%if videosave, close(v2); end
end