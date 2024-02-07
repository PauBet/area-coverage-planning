function [A, fpList] = sidewinder2(startTime, endTime, tobs,...
    inst, sc, target, inroi, olapx, olapy, slewRate, speedUp)
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
%   > inroi:        matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D, in latitudinal 
%                   coordinates [º]
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > speedUp
% 
% Outputs:
%   > A:            cell matrix of the successive instrument observations,
%                   sorted in chronological order.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > fplist:
%   > coverage:
%   > overlap:
%   > makespan:
%   > nfp:
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables
A        = {};
fpList   = [];
if speedUp, resolution = 'lowres';
else, resolution = 'highres'; end
resolution = 'lowres';
amIntercept = false;

% Future work: we need to solve this incompatibility... Anti-meridian and
% visibility:
% Check ROI visible area from spacecraft
vsbroi = visibleroi(inroi, startTime, target, sc);
roi = interppolygon(vsbroi); % interpolate polygon vertices (for improved 
% accuracy)

% Previous anti-meridian intersection check...
ind = find(diff(sort(inroi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind)
    amIntercept = true;
    roi = inroi;
    roi(roi(:, 1) < 0, 1) = roi(roi(:, 1) < 0, 1) + 360;
    [roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2));
end

% Define target area as a polygon
poly1 = polyshape(roi(:, 1), roi(:, 2));

%% Sidewinder heuristics
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
    
    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    [~, ~, ~, ~, tour] = planSidewinderTour2(target, roi, sc, inst, t, ...
        olapx, olapy, fprintc.angle);

    if isempty(tour)
        disp("ROI is not visible from the observer. Define a new observation geometry")
        exit = true;
        continue
    elseif length(tour) == 1
        lon = tour{1}(1); lat = tour{1}(2);
        fpList(end+1) = fprintc;
        tour = {[lon, lat]};
        disp("FOV projection is larger than ROI surface")
        exit = true;
        continue 
    end

    if ~speedUp
        while ~isempty(tour) && ~exit
            % Compute the footprint of each point in the tour successively and
            % subtract the corresponding area from the target polygon
            a = tour{1}; % observation
            tour(1) = []; % delete this observation from the planned tour

            % Check a.m. intercept...
            if a(1) > 180, a(1) = a(1) - 360; end

            % Compute the observation's footprint
            fprintf('Computing %s FOV projection on %s at %s...', inst, ...
                target, cspice_et2utc(t, 'C', 0));
            fprinti = footprint(t, inst, sc, target, resolution, a(1), a(2), 1);
            v2 = fprinti.fovbsight;

            if ~isempty(fprinti.bvertices)
                fprintf('\n')
                % Check a.m. intercept
                if amIntercept
                    aux = fprinti; ind = aux.bvertices(:, 1) < 0;
                    aux.bvertices(ind, 1) = aux.bvertices(ind, 1) + 360;
                    poly2 = polyshape(aux.bvertices);
                else
                    poly2 = polyshape(fprinti.bvertices); % create footprint
                    % polygon
                end

                A{end + 1} = a; % add it in the list of planned observations
                poly1 = subtract(poly1, poly2); % update uncovered area
                if isempty(poly1.Vertices), exit = true; continue; end

                % Save footprint struct
                fpList(end + 1) = fprinti;

                % New time iteration
                if exist("v1", 'var')
                    t = t + tobs + slewDur(v1, v2, slewRate); % future work
                else
                    t = t + tobs;
                end
                v1 = v2;

                % Future work: automated scheduling
                % lastfp = fprinti;
            else
                fprintf(' Surface not reachable\n')
            end
        end
    end
    
    %% Future work: automated scheduling (in-situ)
    % Stop criteria: if the surface of the remaining uncovered roi area is
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

% OK message
fprintf('Sidewinder successfully executed\n')

if ~speedUp
    % Remove first element of fplist (it was just to set the struct fields)
    if ~isempty(fpList)
        fpList(1) = [];
    else
        return;
    end
    % 
    % % Make-span
    % nfp = length(fpList);
    % if length(fpList) > 1
    %     makespan = fpList(end).t - fpList(1).t;
    % else
    %     makespan = tobs;
    % end

    % % ROI coverage percentage
    %[coverage, overlap] = roicoverage(target, polyroi.Vertices, fpList);
else
    A = tour;
    fpList = [];
    % makespan = tobs*length(A);
    % nfp = length(A);
    % coverage = [];
    % overlap = [];
end
end