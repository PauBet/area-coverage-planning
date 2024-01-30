function [roll, pitch, yaw] = scorientation(ax, refframe, t)
% Given the spacecraft axes, compute the yaw, pitch and roll angles that 
% define its orientation.
% Note: this function computes the aforementioned angles by assuming a 
% ZYX Euler rotation.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        [roll, pitch, yaw] = scorientation(axes, refframe, t)
%
% Inputs:
%   > axes: 
%   > refframe:
%   > t:
%
% Outputs:
%   > roll:
%   > pitch:
%   > yaw:

% Compute the tranformation matrix from refframe to the inertial frame
% 'J2000' at time t
rotmat = cspice_pxform(refframe, 'J2000', t);

% Transformation
ax = rotmat*ax;

% Obtain the Euler rotation angles from the rotation matrix. Note: a ZYX
% Euler convention is assumed
[yaw, pitch, roll] = cspice_m2eul(ax, 3, 2, 1);

% Units transformation ([rad] to [deg])
yaw   = yaw*cspice_dpr;
pitch = pitch*cspice_dpr;
roll  = roll*cspice_dpr;
end