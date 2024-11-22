function [A, fpList] = replanningSidewinder(startTime, endTime, tobs,...
    inst, sc, target, inroi, olapx, olapy, slewRate, varargin)
% Target-fixed Boustrophedon decomposition for spacecraft observation
% planning, adapted from [1].
% This function builds a mosaic acquisition for a spacecraft over a
% specified target (ROI) within a given time frame, using a modified
% Boustrophedon decomposition approach. Additionally, it considers the
% instrument cadence, the spacecraft maneuver speed and the target
% visibility.
% This algorithm re-builds the plan after each iteration, updating it to
% align better the allocated tiles to the subsequent observation 
% geometries of the instrument, as opposed to Sidewinder, which is tied to 
% a fixed (the initial) observation geometry.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [A, fpList] = sidewinder(startTime, endTime, tobs,...
%                   inst, sc, target, inroi, olapx, olapy, slewRate, speedUp)
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
%   > inroi:        matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D, in latitudinal 
%                   coordinates [º]
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > slewRate:     rate at which the spacecraft (or instrument platform)
%                   can slew between observations, in [º/s]
% 
% Outputs:
%   > A:            cell matrix of successive instrument observations,
%                   sorted in chronological order.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > fpList:       list of footprint structures detailing the observation
%                   metadata and coverage
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables
A        = {};
fpList   = [];
amIntercept = false;
if nargin == 11 
    resolution = varargin{1};
else
    resolution = 'lowres';
end

% Check ROI visible area from spacecraft
[vsbroi, ~, visibilityFlag] = visibleroi(inroi, startTime, target, sc); % polygon vertices of
% the visible area
if visibilityFlag
    disp("ROI is not visible from the instrument");
    return; 
end 
roi = interppolygon(vsbroi); % interpolate polygon vertices (for improved 
% accuracy)

% Previous anti-meridian intersection check...
ind = find(diff(sort(inroi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind)
    amIntercept = true;
    roi = inroi;
    roi(roi(:, 1) < 0, 1) = roi(roi(:, 1) < 0, 1) + 360; % adjust longitudes
    [roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2)); % sort 
    % coordinates clockwise
end

% [Issue]: We cannot perform visibility and anti-meridian checks
% simultaneously. This means, either we get a sectioned ROI due to
% visibility or anti-meridian, but both may not happen. This is because the
% interpolation function does not work with anti-meridian intercepts.
% [Future work]: Solve this incompatibility

% Define target area as a polygon
poly1 = polyshape(roi(:, 1), roi(:, 2));
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));

%% Replanning Sidewinder heuristics
% The first time iteration is the starting time in the planning horizon
t = startTime;

% Boolean that defines when to stop covering the target area
exit = false;

% Boolean that indicates the first iteration of the algorithm
firstIt = true;

while t <= endTime && ~exit
    
    if firstIt
        % Initial 2D grid layout discretization: the instrument's FOV is going
        % to be projected onto the uncovered area's centroid and the resulting
        % footprint shape is used to set the grid spatial resolution
        gamma(1) = cx; gamma(2) = cy;
        fprintc = footprint(t, inst, sc, target, resolution, ...
            cx, cy, 1); % centroid footprint

        % Initialize struct that saves footprints (sub-structs)
        if t == startTime
            fpList = struct([]);
            for fn = fieldnames(fprintc)'
                fpList(1).(fn{1}) = [];
            end
        end

        % Initialize variables
        gamma0 = [];
        tour = [];
        grid = [];

        fprinti = fprintc;

        firstIt = false;
    end
    
    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    [tour, grid, flag] = rePlanSidewinderTour(target, roi, sc, inst, t, olapx, ...
        olapy, fprinti.angle, cx, cy, tour, grid, gamma, gamma0);

    if flag
        % The previous observation point (gamma) resulted in an
        % unproductive update, leading to a sterile or invalid observation
        % point requiring removal from the tour (which has already been
        % done within rePlanSidewinderTour function). Therefore, we are not
        % going to process the current observation and simply move to the
        % next one
        gamma0 = A{end};
        if isempty(tour)
            break;
        else
            gamma  = tour{1};
        end
        continue; 
    end

    if ~isempty(tour)
        % Process each point of the tour
        [A, tour, fpList, poly1, t] = processObservation(A, tour, ...
            fpList, poly1, t, slewRate, tobs, amIntercept, inst, sc, target, ...
            resolution);

        gamma0 = A{end};
        if isempty(tour)
            break;
        else
            gamma  = tour{1};
        end

        % If polygon is completely covered, break loop. Otherwise, update
        % roi
        fprinti = fpList(end);
        fparea = polyarea(fprinti.bvertices(:, 1), fprinti.bvertices(:, 2));

        % Stop criteria for small uncovered areas (w.r.t. footprint)
        if area(poly1)/fparea < 1e-4, break;
        %if isempty(poly1.Vertices), break;
        else
            [vsbroi, ~, visibilityFlag] = visibleroi(poly1.Vertices, t, ...
                target, sc); % polygon vertices of the visible area
            if visibilityFlag
                disp("ROI is not visible from the instrument");
                break;
            end
            roi = interppolygon(vsbroi);
            poly1.Vertices = roi;
        end
           
    else
        break;
        % For now, the stop criteria is the end of the tour, re-starts are not
        % optimal for the purposes of the scheduling problem
        % [Future work]: automated scheduling (in-situ). Re-starts may be
        % considered, and we will need to define a criteria to prompt those.
    end
end
clear rePlanSidewinderTour;
fprintf('Replanning Sidewinder successfully executed\n')

% Remove first element of fplist (it was just to set the struct fields)
fpList(1) = [];

end