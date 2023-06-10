% Compute neighbours

function [neighbours, poly_next_footprint] = compute_neighbours(target,footprint,poly_footprint,et,steps_low,footprint_func,x_foot, y_foot)

    dx = ones(1,8)*0.6*(max(footprint(1,:))-min(footprint(1,:))); %0.73
    dy = ones(1,8)*0.65*(max(footprint(2,:))-min(footprint(2,:))); %0.85
    ok = zeros(1,8);
    
    %plot(poly_footprint)

    while sum(ok) ~= 8
        neighbours = [target(1) + dx(1)*x_foot(1) , target(2) + dx(1)*x_foot(2);
                      target(1) + dx(2)*x_foot(1) + dy(2)*y_foot(1) , target(2) + dx(2)*x_foot(2) + dy(2)*y_foot(2);
                      target(1) + dy(3)*y_foot(1) , target(2) + dy(3)*y_foot(2);
                      target(1) - dx(4)*x_foot(1) + dy(4)*y_foot(1) , target(2) - dx(4)*x_foot(2) + dy(4)*y_foot(2);
                      target(1) - dx(5)*x_foot(1) , target(2) - dx(5)*x_foot(2);
                      target(1) - dx(6)*x_foot(1) - dy(6)*y_foot(1) , target(2) - dx(6)*x_foot(2) - dy(6)*y_foot(2);
                      target(1) - dy(7)*y_foot(1) , target(2)- dy(7)*y_foot(2);
                      target(1) + dx(8)*x_foot(1) - dy(8)*y_foot(1), target(2) + dx(8)*x_foot(2) - dy(8)*y_foot(2)];
        ok(:) = 1;

        for i = 1:8
            next_fp = footprint_func(neighbours(i,1),neighbours(i,2),et);
            next_footprint = next_fp.bvertices;
            next_footprint = next_footprint.';
            poly_next_footprint{i} = polyshape(next_footprint.');
            poly_out = intersect(poly_footprint,poly_next_footprint{i});
            %plot(poly_next_footprint{i})
            %[cov_lon, cov_lat] = polyclip(next_footprint(1,:),next_footprint(2,:),footprint(1,:),footprint(2,:),'int');  
            if poly_out.NumRegions == 0 %isempty(cov_lon)
                ok(i) = 0;
                dx(i) = dx(i)*0.7;
                dy(i) = dy(i)*0.7;
            end 
        end                   
    end    
end