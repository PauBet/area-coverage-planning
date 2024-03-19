function [fovbounds, boresight, rotmat, visible, varargout] = instpointing(inst, target, sc, t, varargin)
% This function sets the instrument's orientation, provided a target and
% the latitudinal coordinates of the point the instrument should be aiming
% at. It also checks if this point is actually visible from the FOV (could
% be on the dark side of the object as seen from the instrument).
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         01/2023
% 
% Usage:        [fovbounds, boresight, rotmat, visible] = instorient(inst, 
%                   target, sc, t, lon, lat)
%               [fovbounds, boresight, rotmat, visible, lon, lat] = instorient(inst, 
%                   target, sc, t)
%
% Inputs:
%   > inst:       string name of the instrument
%   > target:     string name of the target body
%   > sc:         string name of the spacecraft
%   > t:          observation time, i.e. the minimum time that the 
%                 instrument needs to perform an observation, in seconds
%   > lon:        longitude coordinate of the target body at which the 
%                 instrument boresight is pointing, in [deg]
%   > lat:        latitude coordinate of the target body at which the 
%                 instrument boresight is pointing, in [deg]
%
% Outputs:
%   > fovbounds:  FOV bounds in the body-fixed reference frame centered at 
%                 the spacecraft position at time t
%   > boresight:  FOV boresight in the body-fixed reference frame
%                 centered at the spacecraft position at time t
%   > rotmat:     rotation matrix from instrument frame to target frame
%   > visible:    boolean that determines if the point is visible from the
%                 instrument's FOV
%   > lon:        longitude coordinate of the target body at which the 
%                 instrument boresight is pointing, in [deg]
%   > lat:        latitude coordinate of the target body at which the 
%                 instrument boresight is pointing, in [deg]

% Pre-allocate variables
axis3  = false; % boolean variable that indicates if the spacecraft is 3-axis steerable
if nargin > 4
    lon    = varargin{1};
    lat    = varargin{2};
    axis3  = true;
end
method = 'ELLIPSOID'; % assumption: ray intercept function is going to
% model the target body as a tri-axial ellipsoid
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
abcorr = 'LT'; % one-way light time aberration correction parameter.

% Nargout variables
if nargout>4
    varargout{1} = NaN;
    varargout{2} = NaN;
end

% Retrieve FOV parameters
[shape, instframe, boresight, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % instrument FOV's boundary
    % vectors in the instrument frame
if strcmp(shape, 'CIRCLE') || isequal(shape, 'ELLIPSE')
    error("Circular and ellipsoidal FOV shapes have not been implemented yet")
end
fovbounds = zeros(3, length(bounds));
rotmat = zeros(3);
visible = false;

%% Pointing matrix
if axis3 % 3-axis steerable
    lon = lon*cspice_rpd; lat = lat*cspice_rpd;
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

    % Define a consistent reference vector, e.g., [0; 0; 1] (could be the celestial north)
    reference_vector = [0; 0; 1];

    % Check if boresight is aligned or anti-aligned with reference vector
    if abs(dot(boresight, reference_vector)) > 0.999
        % Adjust reference vector if aligned/anti-aligned to avoid singularity
        reference_vector = [0; 1; 0];
    end

    % Define y-axis using cross product to ensure perpendicularity
    yinst = cross(boresight, reference_vector);
    yinst = normalize(yinst, 'norm'); % Normalize yinst

    % Define x-axis using cross product between boresight and yinst
    xinst = cross(yinst, boresight);
    xinst = normalize(xinst, 'norm'); % Normalize xinst

    % Assign to rotation matrix
    rotmat(:, 1) = xinst;
    rotmat(:, 2) = yinst;
    rotmat(:, 3) = boresight;

else % Pointing constrained by ckernel
    rotmat = cspice_pxform(instframe, targetframe, t);
    [xpoint, ~, ~, found] = cspice_sincpt(method, target, t,...
        targetframe, abcorr, sc, instframe, boresight);
    if found
        [~, lon, lat] = cspice_reclat(xpoint);
        lon = lon*cspice_dpr; lat = lat*cspice_dpr;
    else
        % not visible
        disp("On " + cspice_et2utc(t, 'C', 0) + ", " + inst + " is not pointing at " + target)
        return;
    end
    % Boresight of the instrument must point at the target point
    recpoint = xpoint; % rectangular
    % coordinates of the target point in the body-fixed reference frame
    instpos  = cspice_spkpos(sc, t, targetframe, abcorr, target); % rectangular
    % coordinates of the instrument in the body-fixed reference frame
    v1 = recpoint - instpos; % distance vector to the target point from the
    % instrument in the body-fixed reference frame
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

% Output values
if nargout>4
    varargout{1} = lon;
    varargout{2} = lat;
end
    
end