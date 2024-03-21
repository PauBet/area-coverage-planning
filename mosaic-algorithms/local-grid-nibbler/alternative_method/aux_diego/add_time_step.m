function [ptime, A, fplist, s_area, poly_roi] = add_time_step(time, fpcoverage, poly_roi, poly_target_footprint, A, fplist, target, target_fp, neighbours, neigh_indexes, tobs, inst, target_body, sc, slewRate)

% This function checks if the first target provides enough coverage of the ROI, if
% TRUE it add the target to the output lists A and fplist, and adds a
% timestep (timesteps only added when a valid target is found). The
% variables related to the remaining ROI are updated, and a list of
% possible instants T0 + delta_T is given (each element corresponds to
% one of the 8 possible neighbours)
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date 03/24
%
% Usage: [ptime, A, fplist, s_area, poly_roi] = add_time_step(time, fpcoverage, poly_roi, poly_target_footprint, A, fplist, target, target_fp, ...
%                                                             neighbours, neigh_indexes, tobs, inst, target_body, sc, slewRate)
%
% Inputs:
%   > time:                 Imaging time of the last valid footprint.
%   > fpcoverage:           Coverage of the current target's footprint,
%   > poly_roi:             Polyshape object of the ROI.
%   > poly_target_footprint: Polyshape object of the current target's footprint.  
%   > A:                    cell matrix of successive instrument observations,
%                           sorted in chronological order.
%                           Each observation is defined by the instrument boresight
%                           projection onto the body surface, in latitudinal
%                           coordinates [lon lat], in deg,
%   > fplist:               list of footprint structures detailing the observation
%                           metadata and coverage.
%   > target:               Longitude and latitude of the current target.
%   > target_fp:            Struct containg current footprint features.
%   > neighbours:           array with the valid neighbours coordinates in
%                           longitude and latitude.
%   > neigh_indexes:        list of integers containing the valid neigbours
%                           that can be chosen in the next step.
%   > tobs:                 observation time, i.e. the minimum time that the 
%                           instrument needs to perform an observation, in
%                           seconds.
%   > inst:                 string name of the instrument.
%   > target_body:          string name of the target body.
%   > sc:                   string name of the spacecraft.
%   > slewRate              rate at which the spacecraft (or instrument platform)
%                           can slew between observations, in [º/s].
% Outputs:
%   > ptime:                List of possible imaging times for each valid
%                           neighbour.
%   > A:                    Same as Input with a possible new element
%   > fplist:               Same as Input with a possible new element
%   > s_area:               Scalar containing the remaining uncovered ROI
%                           area.
%   > poly_roi:             Polyshape object containing the remaining
%                           uncovered ROI.


% If the first footprint does not meet the minimum area
ptime = time*ones(1,8);

% If the actual footprint covers part of the real ROI save it and move
% to the next instant               
if fpcoverage >= 0.2
    % Compute the remaining area after substracting the footprint to the
    % ROI
    poly_roi = subtract(poly_roi,poly_target_footprint);
    % Update the area variable
    s_area = area(poly_roi);
    A{end+1} = [target(1), target(2)];
    fplist(end + 1) = target_fp;
    % Possible update times based on neighbours
    for i = neigh_indexes
        ptime(i) = time + tobs + slewDur(target, neighbours(i,:), time, tobs, inst, target_body, sc, slewRate);
    end
else
    s_area = area(poly_roi);
end