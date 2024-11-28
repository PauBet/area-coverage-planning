function [neighbour, neighbours_poly] = compute_neighbours_single(target_body,et,target_orientation,obs,scinst,target_fixed,footprint_func,target,neigh_index)

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
    
    i = neigh_index;
    [spoints, ~, ~, found] = cspice_sincpt('ELLIPSOID', target_body, et, target_fixed, 'NONE', obs,  target_fixed, neigh_directions(:,i)); 
    [~, neighbour(1,1), neighbour(1,2)] = cspice_reclat(spoints);
    neighbour = (180/pi)*neighbour;

    if found == 0
       neighbour = target;
    end    
    
    neigh_fp = footprint_func(neighbour(1,1), neighbour(1,2),et);
    neigh_footprint = neigh_fp.bvertices;
    
    try
        neighbours_poly = polyshape(neigh_footprint);
    catch
        neighbours_poly = 0;
    end    
 
end