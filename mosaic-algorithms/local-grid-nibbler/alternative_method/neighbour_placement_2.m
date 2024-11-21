function [A, fplist, count] = neighbour_placement_2(startTime, tobs, inst, sc, ...
    target_body, roi, olapx, olapy, slewRate)


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

%% Preparing ROI
%Update ROI in case there is an antimeridian cut, transform it to a more convenient way
% Previous anti-meridian intersection check...
ind = find(diff(sort(inroi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind)
    amIntercept = true;
    roi = inroi;
    roi(roi(:, 1) < 0, 1) = roi(roi(:, 1) < 0, 1) + 360; % adjust longitudes
    [roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2)); % sort
    % coordinates clockwise
end

%% Generating tour and choosing first point
%Generating the tour provides with an starting point and a coverage
%path direction
[gamma(1), gamma(2)] = centroid(polyshape(roi(:, 1), roi(:, 2)));
fprintc = footprint(startTime, inst, sc, target_body, 'lowres', ...
    gamma(1), gamma(2), 1); % centroid footprint
tour = planSidewinderTour(target_body, roi, sc, inst, startTime, olapx, olapy, fprintc.angle);
%[tour] = get_tour_first(startTime, inst, sc, ...
%target_body, roi);

[aux, ~, flag] = check_antimeridian_cut(inroi);
if flag
    roi = aux;
end

%The first point in the tour becomes the first target
target = tour{1};

%Create a polyshape object of the ROI that will be updated as it is covered
poly_roi = polyshape(roi);

%% Preliminary steps
%Target reference frame (step related to the use of SPICE)
target_fixed = append('IAU_',target_body);

%Simplify the footprint call for readability
footprint_func = @(x,y,z) footprint_gn(z, inst, sc, target_body, 'lowres', x, y, 1);

%Obtain an example to prepare the fpList struct (allocation)
zerofp = footprint_func(target(1), target(2),startTime);

fplist = struct([]);
for fn = fieldnames(zerofp)'
    fplist(1).(fn{1}) = [];
end


count = 1; % start counting footprint calls (for paper results)

%% Loop parameters
A = []; % Real targets that cover the ROI
error_perc = 0.001; % maximum remaining area acceptable
s_area_zero = area(poly_roi); % Original area
s_area = 10e5; % Remaining area (big value to start)
n = 1; % Number of footprints placed
time = startTime;

%% MOSAIC LOOP

while s_area>s_area_zero*error_perc

    %% Compute the footprint and the neighbours

    if n == 1 % Firstly the coverage direction is defined based on the tour

        % From the target coordinates obtain its footprint and rotation
        % matrix (for neighbours)

        [poly_target_footprint, rotmatrix, target_fp] =  get_poly_rot(target, time, footprint_func);

        count = count + 1;

        % Neighbours coordinates computed before updating time (first iteration all are computed)
        for i = 1:8
            [neighbours(i,:)] = get_neighbour_coordinates(fprintc.angle,target_body,time,rotmatrix,sc,inst,target_fixed,target,i);
        end

        % Compute the coverage of the current footprint over the real ROI
        real_coverage = intersect(poly_roi,poly_target_footprint);
        fpcoverage = area(real_coverage)/area(poly_target_footprint);

        % If the first footprint does not meet the minimum area
        ptime = time*ones(1,8);

        % If the actual footprint covers part of the real ROI save it and move
        % to the next instant
        if fpcoverage >= 0.2
            fprintf('Adding %s FOV projection on %s at %s...\n', inst, ...
                target_body, cspice_et2utc(time, 'C', 0));
            % Compute the remaining area after substracting the footprint to the
            % ROI
            poly_roi = subtract(poly_roi,poly_target_footprint);
            % Update the area variable
            s_area = area(poly_roi);
            A{end+1} = [target(1), target(2)];
            fplist(end + 1) = target_fp;
            % Possible update times based on neighbours
            for i = 1:8
                ptime(i) = time + tobs + slewDur(target, neighbours(i,:), time, tobs, inst, target_body, sc, slewRate);
            end
        end

        % Compute neighbours footprints and its coverage

        for i = 1:8
            [neigh_poly{i}, neigh_rotation_mat(:,:,i), neigh_fp{i}] =  get_poly_rot(neighbours(i,:), ptime(i), footprint_func);
            neigh_coverage(i) = compute_neigh_coverage_single(poly_roi,neigh_poly{i});
            count = count + 1;
        end

        % Based on the coverage path direction defined by the tour and the coverage
        % of each neighbour, we define a list with the 4 neighbours
        % that will be used later to find new targets

        % Compare neighbours with next tour target
        compare_distances = vecnorm((neighbours-tour{2}).');
        compare_distances(1,[2,4,6,8]) = 10^5; % removing diagonal neighbours
        [~, dist_ind] = sort(compare_distances);
        next_target = dist_ind(1);
        neigh_indexes = define_neigh_available(dist_ind, neigh_coverage);

        % Define next target and save the previous one for slew
        % duration
        previous_target = target;
        target = neighbours(next_target, :);

        target_fp = neigh_fp{next_target};
        rotmatrix = neigh_rotation_mat(:,:,next_target);
        poly_target_footprint = neigh_poly{next_target};
        fpcoverage = neigh_coverage(next_target);

        % Update time is required
        if ~isempty(A)
            time = ptime(next_target);
        end

        n = n+1;

    else % In this case only a few neighbours are computed, since by knowing the coverage path, 4 of them are discarded (taboo neighbours
        % )

        for i = neigh_indexes
            [neighbours(i,:)] = get_neighbour_coordinates(fprintc.angle, target_body,time,rotmatrix,sc,inst,target_fixed,target,i);
        end

        % If the first footprint does not meet the minimum area
        ptime = time*ones(1,8);

        % If the actual footprint covers part of the real ROI save it and move
        % to the next instant
        if fpcoverage >= 0.05
            fprintf('Adding %s FOV projection on %s at %s...\n', inst, ...
                target_body, cspice_et2utc(time, 'C', 0));
            % Compute the remaining area after substracting the footprint to the
            % ROI
            poly_roi = subtract(poly_roi,poly_target_footprint);
            % Update the area variable
            s_area = area(poly_roi);
            A{end+1} = [target(1), target(2)];
            fplist(end + 1) = target_fp;
            % Possible update times based on neighbours
            for i = neigh_indexes
                ptime(i) = time + tobs + slewDur(target, neighbours(i,:), time, tobs, inst, target_body, sc, slewRate);
            end
        end

        % Select next target
        if s_area>s_area_zero*error_perc
            [next_target, target_fp, poly_target_footprint, neigh_indexes, rotmatrix, fpcoverage, count] = select_next_target(neigh_indexes,neighbours,poly_roi,ptime,footprint_func,count);
        end

        if next_target == 0
            break
        end

        target = neighbours(next_target,:);
        time = ptime(next_target);

        n = n+1;

    end
end
fprintf('Grid Nibbler successfully executed')
fplist(1) = [];

end