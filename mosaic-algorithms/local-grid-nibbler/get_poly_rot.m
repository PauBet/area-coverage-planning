function [footprint_poly, rotation_mat, target_fp] =  get_poly_rot(target, time, footprint_func, rot_func)
% This function computes the footprint for a given target and time in
% lon-lat, also provides a polyshape version of the footprint (Matlab
% implementation) and the rotation matrix that is necessary to generate the
% coordinates of the neighbours
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date 03/24
%
% Usage: [footprint_poly, rotation_mat, target_fp] =  get_poly_rot(target, time, footprint_func, rot_func)
%
%
% Inputs:
%   > target:       Longitude and latitude of the target on the body
%                   surface aligned with the camera's boresight.
%   > time   :      Imaging time, in TBD seconds past
%                   J2000 epoch
%   > footprint_func:  function handle of the footprint function, so that
%                      only the target and the imaging time is required
%   > rot_func:        function handle of the rotation matrix function, so that
%                      only the target and the imaging time is required
% Outputs:
%   > footprint_poly:    polyshape object with the geometry of the
%                        generated footprint
%   > rotation_mat:      matrix that contains the orientation of the camera
%                        reference frame when pointing to the given target
%                        at the specified time
%   > target_fp:         struct with the footprint features corresponding
%                        to the given imaging time and target
%

    % Compute the footprint given the target and imaging time, other
    % parameters already contained (function handle)
    target_fp = footprint_func(target(1,1), target(1,2),time);
    
    % Compute the corresponding rotation matrix (IAU_TARGET <-> Camera)
    [~, ~, rotation_mat, ~, ~] = rot_func(target(1,1), target(1,2), time);
    
    % Compute the polyshape from the footprint geometry if there is a valid
    % footprint
    try
        footprint_poly = polyshape(target_fp.bvertices);
    catch
        footprint_poly = 0;
    end    

end