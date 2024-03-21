function [neighbour] = get_neighbour_coordinates(angle, target_body,et,target_orientation,obs,scinst,target_fixed,target,neigh_index)
% Computes the neighbour coordinates given the current target, time,
% relevant external conditions, and the desired neighbour index from the 8
% possibilities.
%
% Programmers: Diego Andía (UPC/ESEIAAT)
%              Paula Betriu (UPC/ESEIAAT)
% Date:        08/2023
% Version:      1
% Last update: 03/2024
%
% Usage:       [neighbour] = get_neighbour_coordinates(angle, target_body,et,target_orientation,obs,scinst,target_fixed,target,neigh_index)
%
%
% Inputs: 
%    > angle:                 Rolling angle
%    > target_body:           SPICE identifier of the target body e.g
%                             'EUROPA'
%    > et:                    Current instant, in TDB seconds past
%                             J2000 epoch
%    > target_orientation:    Rotation matrix corresponding to the camera
%                             (from Camera frame to IAU_TARGET)
%                             when pointing the current target at the current instant
%    > obs:                   SPICE identifier of the observer
%    > scinst:                SPICE identifier of the imaging instrumnt
%    > target_fixed:          Body fixed reference frame of the target body
%    > target:                Longitude and latitude of the current target.
%    > neigh_index:           Integer from 1 to 8 that indicates which neighbour to be computed. 
% Outputs:
%    > neighbour:             Longitude and latitude of the obtained neighbour.


% Roll angle around X axis 
rotmat = [cosd(angle)   -sind(angle)  0;
          sind(angle)   cosd(angle)   0;
          0                 0         1];

%Obtains the instrument code
[code_inst, ~] = cspice_bodn2c(scinst); 
%PRECISION REQUIRED FOR THE VALUES
room = 10;
%Obtaining the FOV frame and the bounds vectors.
[~, ~ , ~, bounds] = cspice_getfov(code_inst, room); 
% Define the overlap between neighbours (Future work: Adding it as an external parameter)
overlap = 0.8;

% Obtain a reference distance from the image frame size
dist = 2*bounds(1,1);

% From the reference distance build squared grid of 8 neighbours around the target 
targets = overlap*[-dist, 0, 1/overlap; 
                   -dist, dist, 1/overlap;
                   0, dist, 1/overlap;
                   dist, dist, 1/overlap;
                   dist, 0, 1/overlap; 
                   dist,-dist, 1/overlap;
                   0,-dist, 1/overlap;
                   -dist,-dist, 1/overlap;].';

% Add rolling angle
for i=1:size(targets, 2)
    targets(:, i) = rotmat*targets(:, i);
end

% Convert from the camera reference frame to the IAU target obtaing a vector
% direction
neigh_directions = target_orientation*targets;

% Compute using the corresponding SPICE function the intersection with the
% body surface and express it in degrees. For the specified neighbour
% index.
[spoints, ~, ~, found] = cspice_sincpt('ELLIPSOID', target_body, et, target_fixed, 'NONE', obs,  target_fixed, neigh_directions(:,neigh_index)); 
[~, neighbour(1,1), neighbour(1,2)] = cspice_reclat(spoints);
neighbour(1,:) = (180/pi)*neighbour(1,:);

if found == 0
   neighbour(1,:) = target;
end      

end