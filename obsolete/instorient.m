function [fovbounds, boresight, rotmat, visible] = instorient(inst, target, ...
    lon, lat, sc, t, theta)
% This function sets the instrument's orientation, provided a target and
% the latitudinal coordinates of the point the instrument should be aiming
% at. It also checks if this point is actually visible from the FOV (could
% be on the dark side of the object as seen from the instrument).
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         01/2023
% 
% Usage:        [fovbounds, boresight, rotmat, visible] = instorient(inst, 
%                   target, lon, lat, sc, t, theta)
%
% Inputs:
%   > inst:       string name of the instrument
%   > target:     string name of the target body
%   > lon:        longitude coordinate of the target body at which the 
%                 instrument boresight is pointing, in [deg]
%   > lat:        latitude coordinate of the target body at which the 
%                 instrument boresight is pointing, in [deg]
%   > sc:         string name of the spacecraft
%   > t:          observation time, i.e. the minimum time that the 
%                 instrument needs to perform an observation, in seconds
%   > theta:      rotation angle with respect to the FOV's boresight (DOF),
%                 in deg
%
% Outputs:
%   > fovbounds:  FOV bounds in the body-fixed reference frame centered at 
%                 the spacecraft position at time t
%   > boresight:  FOV boresight in the body-fixed reference frame
%                 centered at the spacecraft position at time t
%   > rotmat:     rotation matrix from instrument frame to target frame
%   > visible:    boolean that determines if the point is visible from the
%                 instrument's FOV

%%
% Pre-allocate variables
boresight = zeros(3, 1);
method = 'ELLIPSOID'; % assumption: ray intercept function is going to 
% model the target body as a tri-axial ellipsoid
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
abcorr = 'LT'; % one-way light time aberration correction parameter.
[~, instframe, boresight, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % instrument FOV's boundary
    % vectors in the instrument frame
fovbounds = zeros(3, length(bounds));
rotmat = zeros(3);
lon = lon*cspice_rpd; lat = lat*cspice_rpd;
visible = false;

%% Instrument orientation calculation
inst2tgt = cspice_pxform(instframe, targetframe, t);
%[~, etto, ~, ~] = cspice_sincpt(method, target, t,...
%    targetframe, 'CN+S', sc, instframe, boresight);
%inst2tgt = cspice_pxfrm2(instframe, targetframe, t, etto);
for i=1:length(bounds) % rotate from instrument to target frame
    bounds(:, i) = inst2tgt*bounds(:, i);
end
boresight = inst2tgt*boresight;

recpoint = cspice_srfrec(cspice_bodn2c(target), lon, lat); % rectangular
% coordinates of the target point in the body-fixed reference frame
instpos  = cspice_spkpos(sc, t, targetframe, abcorr, target); % rectangular
% coordinates of the instrument in the body-fixed reference frame
v2 = recpoint - instpos; % distance vector to the target point from the 
% instrument in the body-fixed reference frame
if dot(v2, recpoint) > 0 % check if the point is visible as seen from the 
    % instrument
    return; 
else
    visible = true;
end
rotAxis = normalize(cross(v2, boresight), 'norm'); % rotation axis over 
% which the instrument pointing has to be rotated, i.e., 
% the cross vector of the body-fixed uz and the final pointing vector (v2)
angle = cspice_vsep(v2, boresight); % phase that has to be rotated from one
% vector to the other
rotmat = cspice_axisar(rotAxis, -angle); % rotation matrix 
rz = cspice_rotate(-theta, 3); % theta rotation (DOF)
rotmat = rotmat*rz; % final rotation matrix
for i=1:length(bounds)
    fovbounds(:,i) = rotmat*bounds(:,i); % instrument FOV's boundary
    % vectors in the target frame (bounds)
end
boresight = rotmat*boresight;

% pointing rotation matrix (and transformation from instrument frame to
% target frame coordinates)
rotmat = rotmat*inst2tgt;

end