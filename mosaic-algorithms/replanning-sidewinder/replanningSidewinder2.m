function [A, fpList, coverage, overlap, makespan, nfp] = replanningSidewinder2(startTime, endTime, tobs,...
    inst, sc, target, inroi, olapx, olapy, speedUp)

% Pre-allocate variables
A        = {};
fpList   = [];
coverage = 0;
overlap  = 0;
makespan = 0;
nfp      = 0;
poly1    = polyshape(inroi(:, 1), inroi(:, 2));
if speedUp, resolution = 'lowres';
else, resolution = 'highres'; end

% Previous anti-meridian intersection check...
ind = find(diff(sort(inroi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind)
    inroi(inroi(:, 1) < 0, 1) = inroi(inroi(:, 1) < 0, 1) + 360;
    [inroi(:, 1), inroi(:, 2)] = sortcw(inroi(:, 1), inroi(:, 2));
end

% Check ROI visible area from spacecraft
vsbroi = visibleroi(inroi, startTime, target, sc);
roi = interppolygon(vsbroi); % interpolate polygon vertices (for improved 
% accuracy)
x = roi(:, 1); y = roi(:, 2);

% Define target area as a polygon
polyroi = polyshape(x, y);

%% Replanning Sidewinder heuristics
% The first time iteration is the starting time in the planning horizon
t = startTime;

% Boolean that defines when to stop covering the target area
exit = false;

% Boolean that indicates the first iteration of the algorithm
firstIt = true;

while ~iszero(roi) && t <= endTime && ~exit
    
    if firstIt
        % Initial 2D grid layout discretization: the instrument's FOV is going
        % to be projected onto the uncovered area's centroid and the resulting
        % footprint shape is used to set the grid spatial resolution
        [gamma(1), gamma(2)] = centroid(polyroi);
        fprintc = footprint(t, inst, sc, target, resolution, ...
            gamma(1), gamma(2)); % centroid footprint

        % Initialize struct that saves footprints (sub-structs)
        if t == startTime
            fpList = struct([]);
            for fn = fieldnames(fprintc)'
                fpList(1).(fn{1}) = [];
            end
        end

        % Initialize variables
        gamma0 = [];
        grid = [];

        % Closest polygon side to the spacecraft's ground track position (this
        % will determine the coverage path in planSidewinderTour)
        %[dir1, dir2] = closestSide(target, sc, t, roi, fprintc.angle);
        dir1 = 'north';
        dir2 = 'east';
        firstIt = false;
        fprinti = fprintc;
    end
    
    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    [tour, grid] = rePlanSidewinderTour(target, roi, sc, inst, t, olapx, ...
        olapy, fprinti.angle, dir1, dir2, grid, gamma, gamma0);

    if ~isempty(tour)
        % Compute the footprint of each point in the tour successively and
        % subtract the corresponding area from the target polygon
        a = tour{1}; % observation
        tour(1) = [];

        % Check a.m. intercept...
        if a(1) > 180, a(1) = a(1) - 360; end

        % Compute the observation's footprint
        fprintf('Computing %s FOV projection on %s at %s...', inst, ...
            target, cspice_et2utc(t, 'C', 0));
        fprinti = footprint(t, inst, sc, target, 'highres', a(1), a(2));

        if ~isempty(fprinti.bvertices)
            fprintf('\n')
            A{end + 1} = a; % add it in the list of planned observations
            poly2 = polyshape(fprinti.bvertices); % create footprint
            % polygon
            poly1 = subtract(poly1, poly2); % update remaining area

            % Save footprint struct
            fpList(end + 1) = fprinti;

            % New time iteration
            t = t + tobs;
            %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work

            currfp = fprinti;
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
            roi   = interpoly1.Vertices;
        end
        gamma0 = a;
        if isempty(tour)
            break;
        else
            gamma  = tour{1};
        end
    else
        break;
        % We should implement a mechanism to re-start the algorithm, in
        % case needed... future work
    end
end
fprintf('Replanning Sidewinder successfully executed\n')

% Remove first element of fplist (it was just to set the struct fields)
fpList(1) = [];

end