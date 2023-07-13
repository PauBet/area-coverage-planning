function neigh_coverage = compute_neighbours_coverage(neighbours, target_orientation, areapoints, footprint_func, et, steps_low)

    for i = 1:8 % For each neighbour
            next_orientation = new_orientation(et,target_orientation.',neighbours(i,:)); 
            next_footprint = footprint_func(et,next_orientation,steps_low);
            neigh_coverage(i) = 0;
            for j = 1:(size(areapoints,2)/2)
                non_zeroes = [];
                non_zeroes = not(areapoints(:,2*j)==0);
                [sub_lon, sub_lat] = polyclip(areapoints(non_zeroes,2*j-1),areapoints(non_zeroes,2*j),next_footprint(1,:),next_footprint(2,:),'int');
                if ~isempty(sub_lon)
                    sub_cov = polyarea(sub_lon{1},sub_lat{1});
                else
                    sub_cov = 0;
                end
                neigh_coverage(i) = neigh_coverage(i) + sub_cov;
            end
    end
    
end