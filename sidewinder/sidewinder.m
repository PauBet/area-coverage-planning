function [A, coverage] = sidewinder(startTime, endTime, tobs, inst, sc, ...
    target, roi, olapx, olapy, videosave)
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
%   > A:            cell matrix of the successive instrument observations,
%                   sorted in chronological order.
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
theta = -20; % temppppppp
[~, targetFrame, ~] = cspice_cnmfrm(target); % body-fixed frame

% Define target area as a polygon
x = roi(:,1); y = roi(:,2);
poly1 = polyshape(x,y);
roiarea = polysurfarea(roi, target); % surface area enclosed by the roi

%% Figure 2.
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex objects

ax2 = mapPlot('vesta-map.png');
fig2 = gcf;
set(gcf,'units','normalized','OuterPosition',[0.4307,0.3144,0.5656,0.6565]);
plot(ax2, polyshape(roi), 'FaceColor', [0.93,0.69,0.13])

% Animation of figure 2
v2 = VideoWriter('topography_map_sidewinder');
v2.FrameRate = 2;
open(v2);
writeVideo(v2,getframe(fig2));

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
    fprintc = footprint(gamma(1), gamma(2), t, inst, sc, target, ...
        theta);   % centroid footprint
    if isempty(fprintc)
        disp("Region of interest not visible from the instrument")
        return; % the footprint is empty because the roi (in this case, 
        % more specifically, the center point of the roi area) is not
        % visible from the instrument FOV. Therefore, the function is
        % exited. This may not be completely correct because the roi area
        % could be large enough so that its centroid is not visible but any
        % other region inside the area is. This is left for future work.
    end

    % Closest polygon side to the spacecraft's ground track position (this
    % will determine the coverage path in planSidewinderTour)
    cside = closestSide(target, sc, t, roi);

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    tour = planSidewinderTour(target, sc, t, roi, fprintc,...
        olapx, olapy, cside, {gamma}, {});

    while ~isempty(tour)
        % Compute the footprint of each point in the tour successively and
        % subtract the corresponding area from the target polygon
        a = tour{1}; % observation
        tour(1) = []; % delete this observation from the planned tour
        A{end + 1} = a; % add it in the list of planned observations
        fprinti = footprint(a(1), a(2), t, inst, sc, target, ...
            theta); % compute the observation's footprint

        if ~isempty(fprinti.bvertices)
            poly2 = polyshape(fprinti.bvertices); % create footprint polygon
            poly1 = subtract(poly1, poly2); % update uncovered area

            %% Plots
            % Footprint plot in Figure 2
            plot(ax2, poly2, 'FaceColor', [0.00,0.45,0.74])
            if length(A) > 1
                plot(ax2, [A{end-1}(1) A{end}(1)], ...
                    [A{end-1}(2) A{end}(2)], 'w-', 'linewidth', 1)
            end
            drawnow
            
            % Figure showing the footprint projection onto the body surface
            footprint3Dprojection(fprinti, videosave)

            % Spacecraft ground track position in Figure 2
            sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
                t, targetFrame, 'NONE', sc);
            [~, sclon, sclat] = cspice_reclat(sctrack);
            plot(sclon*cspice_dpr, sclat*cspice_dpr, 'y^',...
                'MarkerSize', 6)

            % Save Figure 2 frame in the animation
            writeVideo(v2, getframe(fig2));

            % New time iteration
            t = t + tobs;
            %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work

            % Save last footprint
            lastfp = fprinti;
        end
    end

    % Stop criteria: if the surface of the remaining uncovered roi area is
    % smaller than half of the last footprint size, then it is not worth it
    % to start the tour (for the uncovered roi area) again
    if area(polyshape(poly1.Vertices(:, 1), poly1.Vertices(:, 2))) < ...
            0.7*area(polyshape(fprinti.bvertices(:, 1), ...
            fprinti.bvertices(:, 2)))
        exit = true;
    else
        roi = poly1.Vertices;
    end
end

% ROI coverage percentage
coverage = (roiarea - polysurfarea(poly1.Vertices, target))/roiarea;

% End animations
close(v2)
end