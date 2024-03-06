function [neighbours, neighbours_polys] = compute_neighbours_v2(target_body,et,target_orientation,obs,scinst,target_fixed,footprint_func,target,neigh_indexes)

%Obtains the instrument code
[code_inst, ~] = cspice_bodn2c(scinst); 
%PRECISION REQUIRED FOR THE VALUES
room = 10;
%Obtaining the FOV frame and the bounds vectors.
[~, ~ , ~, bounds] = cspice_getfov(code_inst, room); 

overlap = 0.8;

dist = 2*bounds(1,1);

targets = overlap*[-dist, 0, 1/overlap; 
                   -dist, dist, 1/overlap;
                   0, dist, 1/overlap;
                   dist, dist, 1/overlap;
                   dist, 0, 1/overlap; 
                   dist,-dist, 1/overlap;
                   0,-dist, 1/overlap;
                   -dist,-dist, 1/overlap;].';

neigh_directions = target_orientation*targets;

for i = neigh_indexes
    [spoints(:,i), ~, ~, found] = cspice_sincpt('ELLIPSOID', target_body, et, target_fixed, 'NONE', obs,  target_fixed, neigh_directions(:,i)); 
    [~, neighbours(i,1), neighbours(i,2)] = cspice_reclat(spoints(:,i));
    neighbours(i,:) = (180/pi)*neighbours(i,:);
    if found == 0
       neighbours(i,:) = target;
    end    

    neigh_fp = footprint_func(neighbours(i,1), neighbours(i,2),et);
    neigh_footprint = neigh_fp.bvertices;
    try
        neighbours_polys{i} = polyshape(neigh_footprint);
    catch
        neighbours_polys{i} = 0;
    end    
end

 
end