function fpout = fptranslation(fpin, x, y)
% Given a certain footprint, this function performs a translation of the
% footprint to a new center point
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        fpout = fptranslation(fpin, x, y)
%
% Inputs:
%   > target:       struct footprint that wants to be translated. Check
%                   'footprint' function for further details
%   > x:            new longitude point, must be consistent with the
%                   current footprint units ([º], [rad]...)
%   > y:            new latitude point, must be consistent with the
%                   current footprint units ([º], [rad]...)
% 
% Outputs:
%   > fpout:        translated footprints. The struct fields affected are:
%       # olon (boresight observation longitude)
%       # olat (boresight observation latitude)
%       # bvertices (boundary vertices on the lon-lat map)

% Pre-allocate variables...
fpout = fpin;
fpout.olon = x;
fpout.olat = y;

% Translation
dx = x - fpin.olon;
dy = y - fpin.olat;
fpout.bvertices(:, 1) = fpin.bvertices(:, 1) + dx;
fpout.bvertices(:, 2) = fpin.bvertices(:, 2) + dy;

end