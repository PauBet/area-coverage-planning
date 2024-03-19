function [neighbour] = get_neighbour_coordinates(angle, target_body,et,target_orientation,obs,scinst,target_fixed,target,neigh_index)

%
rotmat = [cosd(angle)   -sind(angle)  0;
          sind(angle)   cosd(angle)   0;
          0                 0         1];

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

for i=1:size(targets, 2)
    targets(:, i) = rotmat*targets(:, i);
end

neigh_directions = target_orientation*targets;


[spoints, ~, ~, found] = cspice_sincpt('ELLIPSOID', target_body, et, target_fixed, 'NONE', obs,  target_fixed, neigh_directions(:,neigh_index)); 
[~, neighbour(1,1), neighbour(1,2)] = cspice_reclat(spoints);
neighbour(1,:) = (180/pi)*neighbour(1,:);

if found == 0
   neighbour(1,:) = target;
end      

end