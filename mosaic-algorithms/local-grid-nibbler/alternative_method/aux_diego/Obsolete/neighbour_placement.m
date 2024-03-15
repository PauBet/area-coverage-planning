function [A, fplist, op_count] = neighbour_placement(startTime, tobs, inst, sc, ...
    target_body, roi)    
    
    %Operations counter

    op_count = 0;
    

    % Check if the ROI cuts the antimeridian, the function will always
    % consider the shortest path between two points despite the
    % antimeridian being cutted

    [roi, roi_real, cut_roi] = check_antimeridian_cut(roi);  
 
    %Target reference frame
    
    target_fixed = append('IAU_',target_body);     
    %zerotarget = distantTargetCorner(target_body, sc, startTime, roi).';  
    zerotarget = closestTargetCorner_2(target_body, sc, startTime, roi).'; 

    %Footprint Paula
    
    %footprint_func = @(x,y,z) footprint_2(x, y, z, inst, sc, target_body, 0);
    footprint_func = @(x,y,z) footprint(z, inst, sc, target_body, 'lowres', x, y, 1);

    zerofp = footprint_func(zerotarget(1), zerotarget(2),startTime);
    
    op_count = op_count + 1;

    zerofootprint = zerofp.bvertices.';
    
    fplist = struct([]);
        for fn = fieldnames(zerofp)'
            fplist(1).(fn{1}) = [];
        end

    poly_roi = polyshape(roi);
    remaining_area = poly_roi;

    %% Footprint orientation parameters
    [theta, ~, ~] = foot_axes(zerofootprint);

    %% COMPUTE THE PSEUDO ROI THAT WILL BE USED
    poly_area = polyshape(roi_real);
    [center_area(1), center_area(2)] = centroid(poly_area);
    pseudo_roi = compute_pseudo_area(roi_real,center_area,-theta,cut_roi);
    poly_pseudo_roi = polyshape(pseudo_roi);

    %% START AT THE CLOSER CORNER FROM THE ORIGINAL POINTING
    target = compute_closer(pseudo_roi,zerotarget);
    
    %Compute first footprint
    % zero_orientation = cspice_pxform(inst,target_fixed,startTime);
    % target_orientation = new_orientation(startTime,zero_orientation.',target); 

    %% Loop parameters 
    prev_target = 0; % Neighbour choice of the instant t-1
    prev_prev_target = 0; % Neighbour choice of the instant t-2
    A = []; % Real targets that cover the ROI
    error_perc = 0.001; % maximum remaining area acceptable
    s_area_zero = area(poly_pseudo_roi); % Original area
    s_area = 10e5; % Remaining area (big value to start)
    n = 1; % Number of footprints placed
    % steps_low = 10; % Number of point on each side of the footprints computed during the loop
    time = startTime;

    %% MOSAIC LOOP

    hold on

    while s_area>s_area_zero*error_perc 

        %% Compute the footprint and the neighbours
        %[target_fp, rotmatrix]= footprint_2(target(1), target(2), time, inst, sc, target_body, 0);
        [target_fp, rotmatrix]= footprint(time, inst, sc, target_body, 'lowres', target(1), target(2), 1);
        op_count = op_count + 1;
        target_footprint = target_fp.bvertices.';
        poly_target_footprint = polyshape(target_footprint.');
        [neighbours, neighbours_polys] = compute_neighbours_v2(target_body,time,rotmatrix,sc,inst,target_fixed,footprint_func,target);
        op_count = op_count + 8;
        %plot(poly_target_footprint);

        % Compute the number of subareas after substracting the area by
        % intersection with a footprint the area might end up splitted in
        % different subareas:
             
        % Compute the remaining area after substracting the footprint to the
        % pseudo ROI

        poly_pseudo_roi = subtract(poly_pseudo_roi,poly_target_footprint);
        s_area = area(poly_pseudo_roi);
        
        % Compute the coverage of the current footprint over the real ROI
                
        real_coverage = intersect(poly_roi,poly_target_footprint);
        remaining_area = subtract(remaining_area,poly_target_footprint);

        % If the actual footprint covers part of the real ROI save it and move
        % to the next instant
        if real_coverage.NumRegions ~= 0 %~isempty(real_lon)
            % Footprint plot
            % plot(ax, poly_target_footprint, 'FaceColor', 'b', 'EdgeColor', ...
            %         'b', 'linewidth', 1, 'FaceAlpha', 0.2)
            A{end+1} = [target(1), target(2)];
            fplist(end + 1) = target_fp;
            % if size(A,1) > 1
            %     if abs(A(end-1,1) - A(end,1)) <= 180 % no coverage path -
            %       % a.m. intercept
            %       plot(ax, [A(end-1,1) A(end,1)], [A(end-1,2) ...
            %           A(end,2)], 'w-', 'linewidth', 1)
            %    end
            % end
            % drawnow
            time = time + tobs;
        end
    
        % Compute the coverage of the neighbours 
    
        for i = 1:8         
            try  
                neigh_cov_poly = intersect(poly_pseudo_roi,neighbours_polys{i}); 
                footp_area = area(neighbours_polys{i});
                neigh_coverage(i) = area(neigh_cov_poly)/footp_area;   
            catch
                neigh_coverage(i) = 0;
            end    
         end    

         if max(neigh_coverage)<0.05
                break    
         end    
 
        % Save the two previous choices 
        if n==2
            prev_target = next_target;
        elseif n>2
            prev_prev_target = prev_target;
            prev_target = next_target;
        end 
        
        if n == 26
            stop_here = 1;
        end

        %% 
        if s_area>s_area_zero*error_perc      
            [next_target] = choose_next_target(prev_target,prev_prev_target,neigh_coverage);        
        end
    
        fprintf('Remaining area: %d  PrevPrev target : %d    Prev target: %d    Next target: %d  \n',s_area, prev_prev_target, prev_target, next_target)
        target = neighbours(next_target,:);
        n = n+1;
    end 

fplist(1) = [];

end