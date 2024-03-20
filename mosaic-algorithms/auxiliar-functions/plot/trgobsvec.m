function [obsvec, dist] = trgobsvec(varargin)
% Distance vector between the target point P and the observer at time t
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        obsvec = trgobsvec(srfpoint, t, target, obs)
%               obsvec = trgobsvec(srfpoint, t, target, obs, frame)
%
% Inputs:
%   > srfpoint: target surface point. It can be input either in latitudinal
%               coordinates (in [deg]) or Cartesian coordinates (in [km])
%               with respect to the body-fixed reference frame
%   > t:        time epoch in TDB seconds past J2000 epoch. It can be
%               either a single point in time or a discretized vector of 
%               different time values 
%   > target:   string SPICE name of the target body
%   > obs:      string SPICE name of the observer body
%   > frame:    string SPICE name of the reference frame with respect to
%               which the vector is going to be expressed. If this 
%               variable is not input, the body-fixed reference frame is
%               used by default
%
% Outputs:
%   > obsvec:   observer position vector as seen from the target surface 
%               point in the target body-fixed reference frame, in [km]

% Pre-allocate variables
srfpoint = varargin{1};
t        = varargin{2};
target   = varargin{3};
obs      = varargin{4};
[~, frame, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
abcorr = 'NONE'; % assumption: this function calculates the distance vector
% between the two objects considering the geometric position, no light 
% aberrations are considered (see cspice_spkpos)

% If srfpoint has been input with the latitudinal coordinates, change to
% rectangular
if length(srfpoint) == 2
    srfpoint = srfpoint.*cspice_rpd; % [deg] to [rad]
    srfpoint = cspice_srfrec(cspice_bodn2c(target), srfpoint(1), ...
            srfpoint(2)); % surface point in rectangular coordinates 
            % (body modeled as a tri-axial ellipsoid)
else
    % if srfpoint is input as a 1x3 instead of a 3x1, transpose array
    if size(srfpoint,2) > 1
        srfpoint = srfpoint';
    end
end

% Compute the observer position as seen from the srfpoint
obspos = cspice_spkpos(obs, t, frame, abcorr, target); % compute
% observer position as seen from the target body in the body-fixed ref
% frame
obsvec = obspos - srfpoint; % srfpoint-observer distance vector

% A different reference frame is requested
if nargin > 4
    from = frame;
    to   = varargin{5};
    rotmat = cspice_pxform(from, to, t); % rotation matrix from body-fixed
    % reference frame to the requested one
    if length(t) > 1
        for i=1:length(rotmat)
            obsvec(:, i) = rotmat(:, :, i)*obsvec(:, i);
        end
    else
        obsvec = rotmat*obsvec;
    end
end

% Compute distance
dist = vecnorm(obsvec);
end