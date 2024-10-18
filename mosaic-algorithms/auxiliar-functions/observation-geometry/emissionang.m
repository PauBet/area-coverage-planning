function angle = emissionang(srfpoint, t, target, obs)
% This function returns the phase angle between the target normal to 
% surface and the distance vector to the observer, from the target surface
% point srfpoint, at time t
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        angle = emissionang(srfpoint, t, target, obs)
%
% Inputs:
%   > srfpoint: target surface point. It can be input either in latitudinal
%               coordinates (in [deg]) or Cartesian coordinates (in [km])
%   > t:        time epoch in TDB seconds past J2000 epoch. It can be
%               either a single point in time or a discretized vector of 
%               different time values 
%   > target:   string SPICE name of the target body
%   > obs:      string SPICE name of the observer body
%
% Outputs:
%   > angle:    angle between the normal surface and the distance vector to
%               the observer, in [deg]

% Parameters
method = 'ELLIPSOID'; % assumption: tri-axial ellipsoid modeling of the 
% target body
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE

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
obsvec = trgobsvec(srfpoint, t, target, obs);

% Obtain the outwards surface normal vector
nrmvec = zeros(3,length(t));
for i=1:length(t)
    nrmvec(:,i) = cspice_srfnrm(method, target, t(i), targetframe,...
        srfpoint); % normal to surface 
end

% Angle between the two vectors
angle = cspice_vsep(obsvec, nrmvec);
angle = angle*cspice_dpr; % [rad] to [deg]

end