function [A, fpList] = sidewinder(startTime, endTime, tobs,...
    inst, sc, target, inroi, olapx, olapy, slewRate, speedUp)
% Target-fixed Boustrophedon decomposition for spacecraft observation
% planning, adapted from [1].
% This function builds a mosaic acquisition for a spacecraft over a
% specified target (ROI) within a given time frame, using a modified
% Boustrophedon decomposition approach. Additionally, it considers the
% instrument cadence, the spacecraft maneuver speed and the target
% visibility.
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
%   > speedUp:      boolean flag to run the algorithm in a faster mode
%                   (lower resolution)
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
A        = {}; % observations
fpList   = []; % footprints list
if speedUp, resolution = 'lowres'; % determine resolution based on speedUp flag
else, resolution = 'highres'; end
resolution = 'lowres';
amIntercept = false; % flag for anti-meridian intersection

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

%% Sidewinder algorithm
% The first time iteration is the starting time in the planning horizon
t = startTime;

% Boolean that defines when to stop covering the target area
exit = false;

while ~exit && t < endTime 

    % Initial 2D grid layout discretization: the instrument's FOV is going
    % to be projected onto the uncovered area's centroid and the resulting
    % footprint shape is used to set the grid spatial resolution
    [gamma(1), gamma(2)] = centroid(poly1);
    fprintc = footprint(t, inst, sc, target, resolution, ...
        gamma(1), gamma(2), 1); % centroid footprint

    % Initialize struct that saves footprints (sub-structs)
    if t == startTime
        fpList = struct([]);
        for fn = fieldnames(fprintc)'
            fpList(1).(fn{1}) = [];
        end
    end
    
    % Discretize ROI area (grid) and plan Sidewinder tour based on a
    % Boustrophedon approach
    tour = planSidewinderTour(target, roi, sc, inst, t, olapx, olapy, ...
        fprintc.angle);
    
    % Handle cases where the FOV projection is larger than the ROI area
    if length(tour) == 1
        fpList(end+1) = fprintc;
        disp("FOV projection is larger than ROI surface")
        exit = true;
        continue 
    end
    
    % Process each point of the tour if not in speedUp mode
    if ~speedUp
        while ~isempty(tour) && ~exit
            [A, tour, fpList, poly1, t] = processObservation(A, tour, ...
                fpList, poly1, t, slewRate, tobs, amIntercept, inst, ...
                sc, target, resolution);

            % Polygon completely covered
            fprinti = fpList(end);
            fparea = polyarea(fprinti.bvertices(:, 1), fprinti.bvertices(:, 2));
            
            % Stop criteria for small uncovered areas (w.r.t. footprint)
            if area(poly1)/fparea < 1e-4, exit = true; end
            %if isempty(poly1.Vertices), exit = true; end
        end
    end
    
    % For now, the stop criteria is the end of the tour, re-starts are not
    % optimal for the purposes of the scheduling problem
    % [Future work]: automated scheduling (in-situ). Re-starts may be
    % considered, and we will need to define a criteria to prompt those.
    exit = true;
end

% OK message
fprintf('Sidewinder successfully executed\n')

if ~speedUp
    % Remove first element of fplist (it was just to set the struct fields)
    if ~isempty(fpList)
        fpList(1) = [];
    else
        return;
    end
else
    A = tour;
    fpList = [];
end
end