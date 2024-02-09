function [A, fpList] = frontierRepair2(startTime, endTime, ...
    tobs, inst, sc, target, inroi, olapx, olapy, slewRate, speedUp)

% Pre-allocate variables
A = {}; % List of observations (successive boresight ground track position)
fpList = [];
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
[cx, cy] = centroid(poly1);

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
    fprintc = footprint(t, inst, sc, target, resolution, ...
        gamma(1), gamma(2), 1);   % centroid footprint
    
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

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    [grid, ~, itour, ~, tour, grid_dirx, grid_diry, dir1, dir2] = ...
        planSidewinderTour2(target, roi, sc, inst, t, olapx, olapy, fprintc.angle);
    grid = cellfun(@(c) c', grid, 'UniformOutput', false); % transpose elements

    origin = itour{1};
    while ~isempty(tour)

        a = tour{1}; % current observation point
        old_origin = origin;

        % Update tour
        tour(1) = []; % delete this observation from the planned tour
        itour(1) = [];

        % Gamma and gamma0 in the focal plane
        if ~isempty(tour)
            gamma = tour{1}; % next origin
            origin = itour{1};
        end

        % Check a.m. intercept...
        if a(1) > 180, a(1) = a(1) - 360; end

        % Compute the observation's footprint
        fprintf('Computing %s FOV projection on %s at %s...', inst, ...
            target, cspice_et2utc(t, 'C', 0));
        fprinti = footprint(t, inst, sc, target, resolution, a(1), a(2), 1);

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

            % Intercept footprint
            %polyinter = intersect(poly1, poly2);
            %if area(polyinter)/area(poly2) < 0.2, continue; end
            if isempty(poly1.Vertices), continue; end

            poly1 = subtract(poly1, poly2); % update remaining area
            A{end + 1} = a; % add it in the list of planned observations

            % Save footprint struct
            fpList(end + 1) = fprinti;

            % New time iteration
            if ~isempty(tour)
                p1 = [fprinti.olon, fprinti.olat];
                p2 = [tour{1}(1), tour{1}(2)];
                t = t + tobs + slewDur(p1, p2, t, inst, target, sc, slewRate);
            end
            
            % Update roi
            roi    = interppolygon(poly1.Vertices);
        else
            fprintf(' Surface not reachable\n')
            continue
        end
        
        if isempty(tour)
            break;
        else

            % Update previous grid with the new tile reference (footprint),
            % looking for new potential tiles and/or disposable ones
            [origin, grid, itour, tour] = updateGrid2(roi, itour, ...
                grid, grid_dirx, grid_diry, fprintc.angle, cx, cy, ...
                olapx, olapy, fprinti.angle, dir1, dir2, origin, old_origin, gamma, t, ...
                inst, sc, target);
        end

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
clear updateGrid2;

% Remove first element of fplist (it was just to set the struct fields)
if ~isempty(fpList)
    fpList(1) = [];
end

% End animations
%if videosave, close(v2); end
end