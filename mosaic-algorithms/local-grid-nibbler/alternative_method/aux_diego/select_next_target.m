function [next_target, target_fp, poly_target_footprint, neigh_indexes, rotmatrix, fpcoverage, count] = select_next_target(neigh_indexes,neighbours,poly_roi,ptime,footprint_func, count)
    
    next_target = 0;
    target_fp = [];
    rotmatrix = [];
    poly_target_footprint = [];
    fpcoverage = [];

    min_area = 0.05;
    n = 1;

    for i = neigh_indexes

        %[neighbour, neighbours_poly] = compute_neighbours_single(target_body,ptime,target_orientation,obs,scinst,target_fixed,footprint_func,target,i);
        [neigh_poly{i}, neigh_rotation_mat(:,:,i), neigh_fp{i}] =  get_poly_rot(neighbours(i,:), ptime(i), footprint_func);
        neigh_coverage(i) = compute_neigh_coverage_single(poly_roi,neigh_poly{i});
        
        count = count + 1;
        
        if neigh_coverage(i)>min_area
            next_target = i;
            if n>1
                opposite = get_opposite_neigh(neigh_indexes(1));
                neigh_indexes = [opposite, flip(neigh_indexes(2:end))];
            end 
            
            % Assign outputs
            target_fp = neigh_fp{i};
            rotmatrix = neigh_rotation_mat(:,:,i);
            poly_target_footprint = neigh_poly{i};
            fpcoverage = neigh_coverage(i);

            break
        end    
        n = n + 1;
    end

end