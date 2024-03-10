function [fovbounds, boresight, rotmat, visible] = instpointing(inst, target, ...
    lon, lat, sc, t, theta, axis3)
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

%% Pointing matrix

if axis3
    % Boresight of the instrument must point at the target point
    recpoint = cspice_srfrec(cspice_bodn2c(target), lon, lat); % rectangular
    % coordinates of the target point in the body-fixed reference frame
    instpos  = cspice_spkpos(sc, t, targetframe, abcorr, target); % rectangular
    % coordinates of the instrument in the body-fixed reference frame
    v1 = recpoint - instpos; % distance vector to the target point from the
    % instrument in the body-fixed reference frame
    boresight = normalize(v1, 'norm'); % boresight of the instrument in the
    % body-fixed reference frame

    % The z-axis of the instrument is the boresight
    rotmat(:, 3) = boresight;

    % x- and y-axes must be defined by a "twist" angle
    % when the orientation is given as a parameter, we could define this angle
    % such that the y-axis points upwards, towards the north pole of the body
    xinst = zeros(3, 1);
    yinst = zeros(3, 1);
    yinst(1) = 0;
    yinst(3) = 1;
    xinst(3) = 1;
    yinst(2) = (-boresight(3)/boresight(2))*yinst(3);
    xinst(2) = xinst(3)*(yinst(1)*boresight(3) - yinst(3)*boresight(1))/...
        (boresight(1)*yinst(2) - boresight(2)*yinst(1));
    xinst(1) = -(xinst(2)*boresight(2) + xinst(3)*boresight(3))/boresight(1);
    rotmat(:, 1) = normalize(xinst, 'norm');
    rotmat(:, 2) = normalize(yinst, 'norm');
else
    % Boresight of the instrument must point at the target point
    recpoint = cspice_srfrec(cspice_bodn2c(target), lon, lat); % rectangular
    % coordinates of the target point in the body-fixed reference frame
    instpos  = cspice_spkpos(sc, t, targetframe, abcorr, target); % rectangular
    % coordinates of the instrument in the body-fixed reference frame
    v1 = recpoint - instpos; % distance vector to the target point from the
    % instrument in the body-fixed reference frame
    rotmat = cspice_pxform(instframe, targetframe, t);
end

% Transform coordinates
for i=1:length(bounds)
    fovbounds(:,i) = rotmat*bounds(:,i); % instrument FOV's boundary
    % vectors in the target frame (bounds)
end

% Check if the point is visible as seen from the instrument
if dot(v1, recpoint) > 0 % check if the point is visible as seen from the 
    % instrument
    return; 
else
    visible = true;
end

end