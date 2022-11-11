function map = addfp2map(vertices, map, vlon, vlat)
% Given the polygon vertices, compute and add the values in the lat-lon map 
% that are inside or at the boundary of the former
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        map = footprint2map(vertices, map, vlon, vlat)
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
%                   length(vlat) = size(map,1)
%
% Outputs:
%   > map:      gridded topography map containing the input polygon

% Previous checks...
if isempty(vertices)
    return;
elseif length(vlon) ~= size(map, 2)
    error("length(vlon) must be equal to the number of columns of map")
elseif length(vlat) ~= size(map, 1)
    error("length(vlat) must be equal to the number of rows of map")
end

% For each lat-lon point in the map, if it is inside the footprint polygon,
% add 1
for i=1:size(map,1)
    lat = vlat(i);
    for j=1:size(map,2)
        lon = vlon(j);
        % Check if the [lon, lat] point is within the boundary vertices
        % of the footprint
        % Note: inpolygon also works with disjoint polygons
        if inpolygon(lon, lat, vertices(:,1), vertices(:,2))
            map(i, j) = map(i, j) + 1;
        end
    end
end

end