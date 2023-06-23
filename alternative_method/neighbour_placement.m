function [A, cv, fplist] = neighbour_placement(startTime, tobs, inst, sc, ...
    target_body, roi, ax)    
    
    % Check if the ROI cuts the antimeridian, the function will always
    % consider the shortest path between two points despite the
    % antimeridian being cutted

    max_lon = max(roi(:,1));
    min_lon = min(roi(:,1)); 

    cut_roi = 0;

    if abs(max_lon-min_lon)>200 
        cut_roi = 1;
    end    
    

    %% WORK IN PROGRESS

    if cut_roi == 1
        
        roi_old = roi;

        for i = 1:size(roi,1)
            if roi(i,1)<0
                roi(i,1) = roi(i,1) + 360;
            end    
        end 
        
        roi_real = roi;
        roi = [];

        roi = adapt_roi_cutted(roi_old);
 
    else
        roi_real = roi;
    end    

    
    

    %Target reference frame
    target_fixed = append('IAU_',target_body); 
    
    %FOOTPRINT FUNCTION
    
    cpoint = closestTargetCorner(target_body, sc, startTime, roi);
    zerotarget = cpoint.';    

    %Footprint Paula
    
    footprint_func = @(x,y,z) footprint(x, y, z, inst, sc, target_body, 0);

    zerofp = footprint_func(zerotarget(1), zerotarget(2),startTime);
    zerofootprint = zerofp.bvertices;
    zerofootprint = zerofootprint.';
    
    fplist = struct([]);
        for fn = fieldnames(zerofp)'
            fplist(1).(fn{1}) = [];
        end

    
    poly_roi = polyshape(roi);
    remaining_area = poly_roi;

    %% Footprint orientation parameters
    [theta, x_foot, y_foot] = foot_axes(zerofootprint);

    %% COMPUTE THE PSEUDO ROI THAT WILL BE USED
    poly_area = polyshape(roi_real);
    [center_area(1), center_area(2)] = centroid(poly_area);

    pseudo_roi = compute_pseudo_area(roi_real,center_area,-theta);
    
    if cut_roi == 1

        k = 1;
        pseudo_roi_old = pseudo_roi;
        
        for i = 1:size(pseudo_roi_old,1)
            if pseudo_roi_old(i,1)>180
                pseudo_roi_old(i,1) = pseudo_roi_old(i,1) - 360;
            end    
        end     

        pseudo_roi = [];
        
        pseudo_roi = adapt_roi_cutted(pseudo_roi_old);
        
    end
    
    poly_pseudo_roi = polyshape(pseudo_roi);


    %% START AT THE CLOSER CORNER FROM THE ORIGINAL POINTING

    target = compute_closer(pseudo_roi,zerotarget);
    
    %Compute first footprint
    zero_orientation = cspice_pxform(inst,target_fixed,startTime);
    target_orientation = new_orientation(startTime,zero_orientation.',target); 

    %% Loop parameters 
    prev_target = 0; % Neighbour choice of the instant t-1
    prev_prev_target = 0; % Neighbour choice of the instant t-2
    A = []; % Real targets that cover the ROI
    error_perc = 0.001; % maximum remaining area acceptable
    s_area_zero = area(poly_pseudo_roi); % Original area
    s_area = 10e5; % Remaining area (big value to start)
    n = 1; % Number of footprints placed
    steps_low = 10; % Number of point on each side of the footprints computed during the loop
    time = startTime;

    %% MOSAIC LOOP
    %plot(pseudo_roi(:,1),pseudo_roi(:,2))
    hold on

    while s_area>s_area_zero*error_perc 

        %% Compute the footprint and the neighbours
        [target_fp, rotmatrix]= footprint(target(1), target(2), time, inst, sc, target_body, 0);
        %target_fp = footprint_func(target(1), target(2),time);
        target_footprint = target_fp.bvertices;
        target_footprint = target_footprint.';
        %plot(target_footprint(1,:),target_footprint(2,:))
        poly_target_footprint = polyshape(target_footprint.');
        
        [neighbours, neighbours_polys] = compute_neighbours_v2(target_body,time,rotmatrix,sc,inst,target_fixed,footprint_func);

        % Compute the number of subareas after substracting the area by
        % intersection with a footprint the area might end up splitted in
        % different subareas:
             
        % Compute the remaining area after substracting the footprint to the
        % pseudo ROI

        poly_pseudo_roi = subtract(poly_pseudo_roi,poly_target_footprint);
        s_area = area(poly_pseudo_roi);
        
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


        % Compute the coverage of the current footprint over the real ROI
                
        real_coverage = intersect(poly_roi,poly_target_footprint);
        remaining_area = subtract(remaining_area,poly_target_footprint);

        % If the actual footprint covers part of the real ROI save it and move
        % to the next instant
        if real_coverage.NumRegions ~= 0 %~isempty(real_lon)
            % Footprint plot
            plot(ax, poly_target_footprint, 'FaceColor', 'b', 'EdgeColor', ...
                    'b', 'linewidth', 1, 'FaceAlpha', 0.2)
            A(end+1,:) = [target(1), target(2), time];
            fplist(end + 1) = target_fp;
            if size(A,1) > 1
                if abs(A(end-1,1) - A(end,1)) <= 180 % no coverage path -
                  % a.m. intercept
                  plot(ax, [A(end-1,1) A(end,1)], [A(end-1,2) ...
                      A(end,2)], 'w-', 'linewidth', 1)
               end
            end
            drawnow
            time = time + tobs;
        end
    
        % If the neighbours don't provide coverage finish the loop
        if sum(neigh_coverage) == 0
            fprintf('Loop finished by lack of neighbour coverage')
            break
        end   
 
        % Save the two previous choices 
        if n==2
            prev_target = next_target;
        elseif n>2
            prev_prev_target = prev_target;
            prev_target = next_target;
        end 
    
        %% 
        if s_area>s_area_zero*error_perc      
            [next_target] = choose_next_target(prev_target,prev_prev_target,neigh_coverage);        
        end
    
        fprintf('Remaining area: %d  PrevPrev target : %d    Prev target: %d    Next target: %d  \n',s_area, prev_prev_target, prev_target, next_target)
        target = neighbours(next_target,:);
        %target_orientation = new_orientation(time,target_orientation.',target);  
        n = n+1;
    end 

cv = (area(poly_roi)-area(remaining_area))/area(poly_roi);
fplist(1) = [];

end