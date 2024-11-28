function [A, fplist, count] = neighbour_placement_2(startTime, tobs, inst, sc, ...
    target_body, roi, olapx, olapy, slewRate)
% This function generates a coverage path based on the ROI and
% position/velocity of the spacecraft. From this coverage path the function
% obtains a Boustrophedon like path to cover the ROI, so that it add
% footprints following this path, each new footprint is a "neighbour" of
% the previous one, which means that its correponding pointing direction is
% equal to the one of the previous target plus a default slew.
%
% Programmers: Diego Andía (UPC/ESEIAAT)
%              Paula Betriu (UPC/ESEIAAT)
% Date 03/24
%
% Usage: [A, fplist, count] = neighbour_placement_2(startTime, tobs, inst, sc, ...
%                                         target, roi, olapx, olapy, slewRate)
%
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
%   > roi:        matrix containing the vertices of the ROI polygon. The
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
% Outputs:
%   > A:            cell matrix of successive instrument observations,
%                   sorted in chronological order.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > fpList:       list of footprint structures detailing the observation
%                   metadata and coverage
%


%% Generating tour and choosing first point
% Storing in a different variable the input ROI, the ROI requires a
% particular treatment in order to generate the coverage path, and a
% different one for other operations of this algorithm.
inroi = roi;
% Check ROI visible area from spacecraft
[vsbroi, ~, visibilityFlag] = visibleroi(roi, startTime, target_body, sc); % polygon vertices of
% the visible area
if visibilityFlag
    disp("ROI is not visible from the instrument");
    return;
end
roi = interppolygon(vsbroi); % interpolate polygon vertices (for improved
% accuracy)

ind = find(diff(sort(inroi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind)
    amIntercept = true;
    roi = inroi;
    roi(roi(:, 1) < 0, 1) = roi(roi(:, 1) < 0, 1) + 360; % adjust longitudes
    [roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2)); % sort
    % coordinates clockwise
end

% Generate a first footprint to obtain a reference regarding angles
[gamma(1), gamma(2)] = centroid(polyshape(roi(:, 1), roi(:, 2)));
fprintc = footprint(startTime, inst, sc, target_body, 'lowres', ...
    gamma(1), gamma(2), 1); % centroid footprint

% Generate the first tour or coverage path that will be used as a
% reference to define the Boustrophedon-like path that the algorithm
% will follow to add neighbours
tour = planSidewinderTour(target_body, roi, sc, inst, startTime, olapx, olapy);

%% Preparing ROI
% Check if the ROI cuts the antimeridian and adapts the ROI
% coordinates. This is only necessary for this implementation of the
% algorithm since it relies on certain Matlab functions (polyshape objects) for polygons
% and it considers the body surface as a longitude-latitude map from
% -180 to +180.
[aux, flag] = check_antimeridian_cut(inroi);
if flag
    roi = aux;
end

it = 0;
%A  = {};
%while it <= length(tour) && isempty(A)
    it = it + 1;

    % The first objective is the first point in the tour
    target = tour{it};

    % Create a polyshape object of the ROI that will be updated as it is
    % covered.
    poly_roi = polyshape(roi);

    %% Pre-allocation of variables and definition of useful functions
    %Target reference frame (step related to the use of SPICE)
    target_fixed = append('IAU_',target_body);

    %Simplify the footprint call for readability(some inputs will be always the same e.g inst, sc, target_body, resolution)
    footprint_func = @(x,y,z) footprint(z, inst, sc, target_body, 'lowres', x, y, 0);

    rot_func = @(x,y,z) instpointing(inst, target_body, sc, z, x, y);

    %Obtain an example to prepare the fpList struct (allocation)
    zerofp = footprint_func(target(1), target(2),startTime);

    %Preallocate the fplist that will be returned using a reference fp
    %struct (OUTPUT)
    fplist = struct([]);
    for fn = fieldnames(zerofp)'
        fplist(1).(fn{1}) = [];
    end
    % Real targets that cover the ROI (OUTPUT)
    A = [];

    count = 1; % start counting footprint calls (for paper results)

    %% Loop parameters
    error_perc = 0.001; % maximum remaining area acceptable
    s_area_zero = area(poly_roi); % Original area
    s_area = 10e5; % Remaining area (big value to start)

    n = 1; % Step counter

    % Time variable that will be updated after adding a new footprint
    time = startTime;


    %% Definition of the order of addition of neighbours

    % For the first target, obtains its corresponding footprint polyshape,
    % the rotation matrix the camera reference system, and its
    % corresponding fp struct
    [poly_target_footprint, rotmatrix, target_fp] =  get_poly_rot(target, ...
        time, footprint_func, rot_func);

    count = count + 1;

    % Compute the coordinatees in lon-lat of the 8 neighbours of the first
    % target. The coordinates are computed at the imaging time of the
    % target itself.
    for i = 1:8
        [neighbours(i,:)] = get_neighbour_coordinates(fprintc.angle,target_body,time,rotmatrix,sc,inst,target_fixed,target,i);
    end

    % Compute the coverage of the footprint corresponding to the first taget over the real ROI
    fpcoverage = get_fp_coverage(poly_roi, poly_target_footprint);

    % Checks if the first target provides enough coverage of the ROI, if
    % TRUE it add the target to the output lists A and fplist, and adds a
    % timestep (timesteps only added when a valid target is found). The
    % variables related to the remaining ROI are updated, and a list of
    % possible instants T0 + delta_T is given (each element corresponds to
    % one of the 8 possible neighbours)
    [ptime, A, fplist, s_area, poly_roi] = add_time_step(time, fpcoverage, poly_roi, poly_target_footprint, A, fplist, target, target_fp, neighbours, [1,2,3,4,5,6,7,8], tobs, inst, target_body, sc, slewRate);
%end

% Compute neighbours footprints and its coverage using the coordinates
% previously obtained
for i = 1:8
    [neigh_poly{i}, neigh_rotation_mat(:,:,i), neigh_fp{i}] =  get_poly_rot(neighbours(i,:), ptime(i), footprint_func, rot_func);
    neigh_coverage(i) = get_fp_coverage(poly_roi,neigh_poly{i});
    count = count + 1;
end

% Find which of the neighbours is closer to the next tour (coverage path previously computed)
compare_distances = vecnorm((neighbours-tour{it+1}).');
% removing diagonal neighbours, this is because the coverage path
% should be followed adding neighbours that are not in diagonal
% positions since they provide less overlap with the current target and
% the likeliness of having gaps between footprints increases
compare_distances(1,[2,4,6,8]) = 10^5;
% The non-diagonal neighbour closer to the next tour direction will be
% chosen as the next target
[~, dist_ind] = sort(compare_distances);
next_target = dist_ind(1); % Number of the neighbour chosen as next target
%% [Paula]: This is not correct in some cases... What if you have a case
% like this:
% X O X
% O O O
% where X corresponds to non-covering and O to covering footprints.
% According to this rationale, the algorithm will choose the
% element in row 2 column 2, however, it should choose either (2, 1) or 
% (2, 3). This needs to be corrected!

% Knowing which neighbour aligns better with the coverage path, and the
% coverage of the rest of neighbours it is possible to identify in
% which order should be check neighbours when covering the ROI, this
% helps to avoid making polygon operations for neighbours that are
% always invalid given the coverage path.
neigh_indexes = define_neigh_available(dist_ind, neigh_coverage);

% Define next target and save the previous one for slew
% duration, assign all the variables related to the next target to the
% chosen neighbour
target = neighbours(next_target, :);
target_fp = neigh_fp{next_target};
rotmatrix = neigh_rotation_mat(:,:,next_target);
poly_target_footprint = neigh_poly{next_target};
fpcoverage = neigh_coverage(next_target);


%% MOSAIC LOOP
% Add footprints until the ROI has been sufficiently covered
flag = 0;
while s_area>s_area_zero*error_perc

    % For each of the valid neighbours saved in neigh_indexes
    % compute the corresponding lon-lat coordinates at the current
    % target's imaging time (To ensure overlap)
    for i = neigh_indexes
        [neighbours(i,:)] = get_neighbour_coordinates(fprintc.angle, ...
            target_body,ptime(next_target),rotmatrix,sc,inst,...
            target_fixed,target,i);
    end

    % Update the time variable with the last valid footprint added,
    % since footprints that do not provide coverage do not add to
    % the makespan, and are only used as a bridge between valid
    % footprints when creating the mosaic.
    time = ptime(next_target);

    % Check if the new target provides enough coverage, if TRUE
    % update possible imaging times of the neighbours, add new
    % target to the output variables, and update remaining ROI
    % related variable
    
    [ptime, A, fplist, s_area, poly_roi] = add_time_step(time, ...
        fpcoverage, poly_roi, poly_target_footprint, A, fplist, ...
        target, target_fp, neighbours, neigh_indexes, tobs, inst, ...
        target_body, sc, slewRate);
    
    % Stop criteria for small uncovered areas (w.r.t. footprint)
    fp = fplist(end);
    fparea = area(polyshape(fp.bvertices));
    if s_area/fparea < 1e-4, break; end

    % Check roi visibility
    [vsbroi, ~, visibilityFlag] = visibleroi(poly_roi.Vertices, ...
        time, target_body, sc);
    if visibilityFlag
        disp("ROI no longer reachable");
        break;
    else
        poly_roi.Vertices = interppolygon(vsbroi);
    end

    % Select next target
    if s_area>s_area_zero*error_perc
        [next_target, target_fp, poly_target_footprint, neigh_indexes, ...
            rotmatrix, fpcoverage, count] = select_next_target(neigh_indexes,neighbours,poly_roi,ptime,footprint_func,rot_func,count);
    end

    % If no valid neighbour is found (Dead-end)
    if next_target == 0
        break
    end

    % Assign the chosen neighbour coordinates to the new target
    % variable
    target = neighbours(next_target,:);
end

% The first element of the fplist is empty so it is removed
fplist(1) = [];

% OK message
fprintf('Grid Nibbler successfully executed\n')

end