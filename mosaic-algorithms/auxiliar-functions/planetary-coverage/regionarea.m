function [A, coverage] = regionarea(body, lon, lat)
% Provided a specific body target and the vertices that bound a region on
% its surface, this function calculates the area of the region.
% Body is modeled as an oblate spheroid.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         12/2023
% Revision:     1
%
% Usage:        [A, coverage] = regionarea(body, lon, lat)
%
% Inputs:
%   > body:     SPICE ID (int) or name (string) of the body
%   > lon:      Longitude values of the vertices that bound the region of
%               interest, in [deg]
%   > lat:      Latitude values of the vertices that bound the region of
%               interest, in [deg]
%
% Output:
%   > A:        Area of the surface, in [km2]. If there are more than one
%               regions, A is an array
%   > coverage: Percentage of surface this area covers with respect to the
%               entire body surface, in [%]. If there are more than one
%               regions, coverage is an array

% Pre-allocate variables
A = 0;

% Previous check...
if isempty(lat) || isempty(lon)
    return;
end

% Retrieve body radii from SPICE
radii = cspice_bodvrd(body,'RADII',3);

% Body modelization as an oblate spheroid, i.e. radii(1) == radii(2)
customPlanet = referenceEllipsoid;
customPlanet.LengthUnit = 'kilometer';
customPlanet.SemimajorAxis = radii(1);
customPlanet.SemiminorAxis = radii(3);

% Surface area
A = areaint(lat, lon, customPlanet, 'degrees');
A = sum(A); % in case polygon has more than one region

% Coverage
GA = customPlanet.SurfaceArea; % total surface area
coverage = A*1e2/GA; % [%]
end