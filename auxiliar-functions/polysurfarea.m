function area = polysurfarea(varargin)
% Given a set of boundary polygon vertices, in latitudinal coordinates, 
% compute the body surface area that it encloses
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        area = fparea(vertices, target)
%               area = fparea(vertices, target, map, vlon, vlat)
%
% Inputs:
%   > vertices: matrix (N, 2) containing the N boundary vertices of the
%               polygon, in latitudinal coordinates, in deg. vertices may
%               contain more than one polygon, split because the actual
%               footprint intersects the anti-meridian of the body surface.
%               In this case, the polygons are separated in 'vertices' by
%               [NaN, NaN]
%               Example: area1 = [350 30; 360 30; 360  0; 350 0]; 
%                        area2 = [ 0 30; 10 30; 10  0; 0  0];
%                         area = [area1; [NaN NaN]; area2];
%                       figure
%                       plot(polyshape(area(:,1), area(:,2)))
%   > target:   SPICE string name of the target body
%
%   Note: The following variables may not be input. In that case, vlon and
%   vlat are two arrays covering spans of [-180, 180]º and [90, -90]º and 1
%   deg of resolution. 'map' is a matrix length(vlat)xlength(vlon) of zeros
%   > map:      int matrix of the surface map of a body. Each
%               element in the map corresponds to a set of [lat, lon]
%               coordinates:
%                   longitude (-) ------> (+)
%               latitude (+)  map(1,1) map(1,2) ...
%                         ¦   map(2,1)
%                         ¦      ⁝
%                         ¦
%                         ∨
%                        (-)
%               Each value in the matrix accounts for the number of passes
%               an instrument has performed over a certain point of the 
%               body surface. Thus, map(p,q) = 1 means that the point
%               P(vlon(q),vlat(p)) has been framed once, map(p,q) 
%               = 2 twice, and so on.
%   > vlon:     array of discretized longitude values, in deg. 
%                   length(vlon) = size(map, 2)
%   > vlat:     array of discretized latitude values, in deg.
%                   length(vlat) = size(map, 1)
%
% Outputs:
%   > fp:       surface area value, in [km2]

% Variables allocation
if nargin < 2
    error("Incorrect number of inputs. The function needs the " + ...
        "polygon vertices and target SPICE name at least")
else
    vertices = varargin{1}; target = varargin{2};
    if nargin == 5
        map = varargin{3}; vlon = varargin{4}; vlat = varargin{5};
    elseif nargin == 2
        vlon = -180:1:180; vlat = 90:-1:-90;
        map = zeros(length(vlat), length(vlon));
    else
        error("Incorrect number of inputs")
    end
end
area = 0;

% Previous checks...
if isempty(vertices)
    return;
elseif length(vlon) ~= size(map, 2)
    error("length(vlon) must be equal to the number of columns of map")
elseif length(vlat) ~= size(map, 1)
    error("length(vlat) must be equal to the number of rows of map")
end

% In case 'vertices' is not empty but it was not added in the surface
% map...
if ~all(map, 'all'), map = addfp2map(vertices, map, vlon, vlat); end

% Parameters
steplon = (vlon(2) - vlon(1))*cspice_rpd; % vlon is uniformly gridded
steplat = (vlat(2) - vlat(1))*cspice_rpd; % vlat is uniformly gridded

% Compute surface area by iterating and adding each surface element
for i=1:size(map,1)
    lat = vlat(i)*cspice_rpd;
    for j=1:size(map,2)
        lon = vlon(j)*cspice_rpd;
        if map(i,j)
            r2 = vecnorm(cspice_srfrec(cspice_bodn2c(target),...
                lon, lat))^2;
            area = area + abs(r2*steplon*steplat*cos(lat));
        end
    end
end

end