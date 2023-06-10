function [A] = neighbour_placement(startTime, tobs, inst, sc, ...
    target_body, roi, ax)
    
    %Target reference frame
    target_fixed = append('IAU_',target_body); 
    
    %FOOTPRINT FUNCTION
    steps_zero = 90; % Number of points for the footprint function

    footprint_func = @(x,y,z) footprint_coverage_final(x, sc, inst, y, target_body, target_fixed,z); 

    %ZERO CONDITIONS
    [zerofootprint,~, ] = footprint_func(startTime,0,steps_zero);
    zerotarget = [mean(zerofootprint(1,:));
             mean(zerofootprint(2,:))];

    %Footprint Paula
    
    footprint_func = @(x,y,z) footprint(x, y, z, inst, sc, target_body, 0);

    zerofp = footprint_func(zerotarget(1), zerotarget(2),startTime);
    zerofootprint = zerofp.bvertices;
    zerofootprint = zerofootprint.';

    roi_real = roi;
    poly_roi = polyshape(roi_real);

    %% Footprint orientation parameters
    [theta, x_foot, y_foot] = foot_axes(zerofootprint,steps_zero);

    %% COMPUTE THE PSEUDO ROI THAT WILL BE USED
    poly_area = polyshape(roi_real);
    [center_area(1), center_area(2)] = centroid(poly_area);

    pseudo_roi = compute_pseudo_area(roi_real,center_area,theta);
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
    s_area_zero = polyarea(pseudo_roi(:,1),pseudo_roi(:,2)); % Original area
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
            neigh_cov_poly = intersect(poly_pseudo_roi,neighbours_polys{i}); 
            neigh_coverage(i) = area(neigh_cov_poly);   
        end    


        % Compute the coverage of the current footprint over the real ROI
                
        real_coverage = intersect(poly_roi,poly_target_footprint);

        % If the actual footprint covers part of the real ROI save it and move
        % to the next instant
        if real_coverage.NumRegions ~= 0 %~isempty(real_lon)
            % Footprint plot
            plot(ax, poly_target_footprint, 'FaceColor', 'b', 'EdgeColor', ...
                    'b', 'linewidth', 1, 'FaceAlpha', 0.2)
            A(end+1,:) = [target(1), target(2), time];
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
    
        if s_area>s_area_zero*error_perc      
            [next_target] = choose_next_target(prev_target,prev_prev_target,neigh_coverage);        
        end
    
        fprintf('Remaining area: %d  PrevPrev target : %d    Prev target: %d    Next target: %d  \n',s_area, prev_prev_target, prev_target, next_target)
        target = neighbours(next_target,:);
        %target_orientation = new_orientation(time,target_orientation.',target);  
        n = n+1;
    end 

end