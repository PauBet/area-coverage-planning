function [A, fpList] = frontierRepair(startTime, endTime, ...
    tobs, inst, sc, target, inroi, olapx, olapy, slewRate, varargin)
% This function adjusts the observation grid and planning tour in response 
% to new observations. It takes into account changes in observation 
% geometry, incorporating new observation points, and removes points that 
% are no longer necessary for covering the region of interest (ROI).
% A Boustrophedon pattern is used to ensure efficient coverage. It has been
% adapted from [1]
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [A, fpList] = frontierRepair(startTime, endTime, ...
%                   tobs, inst, sc, target, inroi, olapx, olapy, slewRate, ...
%                   resolution)
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
%   > resolution:   string 'lowres' or 'highres' that determines the
%                   footprint resolution calculation. It's set to 'lowres'
%                   by default
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
A = {}; % List of observations (successive boresight ground track position)
fpList = [];
amIntercept = false;
if nargin == 11 
    resolution = varargin{1};
else
    resolution = 'lowres';
end

% % Check ROI visible area from spacecraft
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
[cx, cy] = centroid(poly1);

%% Frontier Repair algorithm
% The first time iteration is the starting time in the planning horizon
t = startTime;

% Boolean that defines when to stop covering the target area
exit = false;

while ~exit

    % Initial 2D grid layout discretization: the instrument's FOV is going
    % to be projected onto the uncovered area's centroid and the resulting
    % footprint shape is used to set the grid spatial resolution
    [gamma(1), gamma(2)] = centroid(polyshape(roi(:,1),roi(:,2)));
    fprintc = footprint(t, inst, sc, target, resolution, ...
        gamma(1), gamma(2), 0);   % centroid footprint;
    
    % Initialize struct that saves footprints (sub-structs)
    if t == startTime
        fpList = struct([]);
        for fn = fieldnames(fprintc)'
            fpList(1).(fn{1}) = [];
        end
    end

    % Check roi visibility
    [vsbroi, ~, visibilityFlag] = visibleroi(roi, t, target, sc);
    if visibilityFlag
        disp("ROI no longer reachable");
        break;
    else
        roi = interppolygon(vsbroi);
    end

    % Discretize ROI area (grid) and plan Sidewinder tour based on a
    % Boustrophedon approach
    [tour, grid, itour, grid_dirx, grid_diry, dir1, dir2] = ...
        planSidewinderTour(target, roi, sc, inst, t, olapx, olapy);
    grid = cellfun(@(c) c', grid, 'UniformOutput', false); % transpose elements

    % Handle cases where the FOV projection is larger than the ROI area
    if length(tour) < 1
        A{end + 1} = gamma;
        fpList(end+1) = fprintc;
        exit = true;
        continue
    end

    seed = itour{1};
    while ~isempty(tour) && t < endTime

        % Update origin and tour
        old_seed = seed;
        itour(1) = [];

        % Process each point of the tour
        [A, tour, fpList, poly1, t] = processObservation(A, tour, ...
            fpList, poly1, t, slewRate, tobs, amIntercept, inst, sc, target, ...
            resolution);

        % If polygon is completely covered, break loop
        if isempty(poly1.Vertices), break; end

        % Update roi
        roi = poly1.Vertices;

        % Check roi visibility
        [vsbroi, ~, visibilityFlag] = visibleroi(roi, t, target, sc);
        if visibilityFlag 
            disp("ROI no longer reachable");
            break;
        else
            roi = interppolygon(vsbroi);
        end
        
        if isempty(tour)
            break;
        else
            gamma = tour{1}; % next observation point
            seed = itour{1}; % next seed in the image plane

            % Update previous grid with the new tile reference (footprint),
            % looking for new potential tiles and/or disposable ones
            [seed, grid, itour, tour] = updateGrid(roi, itour, ...
                grid, grid_dirx, grid_diry, cx, cy, ...
                olapx, olapy, dir1, dir2, seed, old_seed, gamma, t, ...
                inst, sc, target);
        end

    end

    % For now, the stop criteria is the end of the tour, re-starts are not
    % optimal for the purposes of the scheduling problem
    % [Future work]: automated scheduling (in-situ). Re-starts may be
    % considered, and we will need to define a criteria to prompt those.
    exit = true;
end
clear updateGrid;
clear checkTaboo;

% OK message
fprintf('Online Frontier successfully executed\n')

% Remove first element of fplist (it was just to set the struct fields)
if ~isempty(fpList)
    fpList(1) = [];
end

end