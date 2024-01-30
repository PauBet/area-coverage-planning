function roi = interppolygon(roi0)
% This function interpolates a polygon defined by longitude and latitude
% points. The interpolation distance is defined according to the minimum
% Euclidean distance between the points that enclose the polygon
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         1/2024
%
% Usage:        roi = interppolygon(roi0)
%
% Input:
%   > roi0:     A Nx2 matrix where each row represents a point in 2D space, 
%               typically [longitude, latitude]
%
% Output:
%   > roi:      Updated roi with interpolated coordinates

% Close polygon
% Check if the first and the last point of the polygon are the same.
% If not, add the first point to the end of the list to close the polygon.
if any(roi0(1, :) ~= roi0(end, :)), roi0(end+1, :) = roi0(1, :); end

% Definition of maximum allowable distance
% Initialize variables to store longitude and latitude from the polygon.
lon = roi0(:, 1);
lat = roi0(:, 2);
epsilon = inf; % initialize epsilon to infinity. This will be used to 
% find the minimum distance between points.

% Loop through each pair of points to find the minimum non-zero distance
for i=1:length(lon)-1
    p1 = [lon(i), lat(i)]; % current point
    p2 = [lon(i+1), lat(i+1)]; % next point
    dist = vecnorm(p2 - p1); % calculate Euclidean distance between points
    if dist == 0, continue; end % skip if points are identical
    if dist < epsilon; epsilon = dist; end % update minimum distance
end

% Perform interpolation of the latitude and longitude based on the 
% minimum distance found
[newlat, newlon] = interpm(lat, lon, ceil(epsilon/2));

% Update the roi array with the interpolated latitude and longitude values
roi(:, 2) = newlat; roi(:, 1) = newlon;

% Sort coordinates in clockwise order
%[roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2));
end