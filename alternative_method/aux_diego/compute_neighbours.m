% Compute neighbours

function neighbours = compute_neighbours(target,footprint,et,steps_low,target_orientation,footprint_func,x_foot, y_foot)

    dx = 0.95*(max(footprint(1,:))-min(footprint(1,:))); %0.73
    dy = 0.95*(max(footprint(2,:))-min(footprint(2,:))); %0.85
    ok = 0;
    while ok == 0
        neighbours = [target(1) + dx*x_foot(1) , target(2) + dx*x_foot(2);
                      target(1) + dx*x_foot(1) + dy*y_foot(1) , target(2) + dx*x_foot(2) + dy*y_foot(2);
                      target(1) + dy*y_foot(1) , target(2) + dy*y_foot(2);
                      target(1) - dx*x_foot(1) + dy*y_foot(1) , target(2) - dx*x_foot(2) + dy*y_foot(2);
                      target(1) - dx*x_foot(1) , target(2) - dx*x_foot(2);
                      target(1) - dx*x_foot(1) - dy*y_foot(1) , target(2) - dx*x_foot(2) - dy*y_foot(2);
                      target(1) - dy*y_foot(1) , target(2)- dy*y_foot(2);
                      target(1) + dx*x_foot(1) - dy*y_foot(1), target(2) + dx*x_foot(2) - dy*y_foot(2)];
        ok = 1;

        for i = 1:8
            next_orientation = new_orientation(et,target_orientation.',neighbours(i,:)); 
            next_footprint = footprint_func(et,next_orientation,steps_low);
            [cov_lon, cov_lat] = polyclip(next_footprint(1,:),next_footprint(2,:),footprint(1,:),footprint(2,:),'int');
            if isempty(cov_lon)
                ok = 0;
            end 
        end    
        
        if ok == 0
                dx = dx*0.9;
                dy = dy*0.95;
        end       
    end    
end