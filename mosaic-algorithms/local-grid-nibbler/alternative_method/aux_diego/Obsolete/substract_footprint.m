function [areapoints, s_area] = substract_footprint(areapoints, target_footprint)
    
    s_area = 0;
    n_sub_areas = size(areapoints,2)/2;
    areapoints_old = areapoints;
    areapoints = zeros(5000,1);
    k = 1;

    for i = 1:n_sub_areas % For each subarea
        % list with the elements of the current subareas that are not empty
        correct_vals_old = not(areapoints_old(:,2*i)==0);
        % substract the footprint coverage to the current subarea 
        [areapoints_lon, areapoints_lat] = polyclip(areapoints_old(correct_vals_old,2*i-1),areapoints_old(correct_vals_old,2*i),target_footprint(1,:),target_footprint(2,:),'difference');
        
        for j = 1:size(areapoints_lon,2)
            size_sub = size(areapoints_lon{j},1); % Compute the size of the resulting subareas
            areapoints(1:size_sub,2*k-1) = areapoints_lon{j}; % Add the subareas to a global subarea array
            areapoints(1:size_sub,2*k) = areapoints_lat{j};
            correct_vals = not(areapoints(:,2*k)==0); % Check the non-empty values of the footprint
            s_area_add = polyarea(areapoints(correct_vals,2*k-1),areapoints(correct_vals,2*k)); % Compute the area of the subarea
            s_area = s_area + s_area_add; % Compute the remaining area to cover
            k = k + 1;
        end    
    end 
end