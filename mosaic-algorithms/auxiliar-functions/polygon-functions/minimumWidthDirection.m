function [thetamin, minwidth, height, axes] = minimumWidthDirection(x, y)
% This function computes the orientation (angle) at which the width of the 
% polygon is minimized. It also calculates the minimum width, the height 
% (assuming the minimum width as the base), and the direction vector (axes) 
% of the minimum width
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        [thetamin, minwidth, height, axes] = minimumWidthDirection(x, y)
%
% Inputs:
%   > x:            array of x-coordinates of the polygon vertices
%   > y:            array of y-coordinates of the polygon vertices
% 
% Outputs:
%   > thetamin:     angle at which polygon's width is minimized, in [deg]
%   > minwidth:     minimum width of the polygon
%   > height:       height of the polygon, at the orientation specified by
%                   thetamin and assuming minwidth as the base
%   > axes:         unit vector representing the direction of the polygon
%                   axes (width, height)

% Check if the polygon is divided in two (a.m. intersection)...
ind = find(isnan(x));
if ~isempty(ind)
    x(1:ind) = x(1:ind) + 360;
    x(ind) = [];
    y(ind) = [];
end

% Find centroid
[cx, cy] = centroid(polyshape(x, y));

% Sort the vertices in clockwise direction
[vertices(:,1), vertices(:,2)] = sortcw(x, y);
%vertices(end+1, :) = vertices(1, :); % close polygon

% Minimum width direction
npoints = 361;
angle = linspace(0, 360, npoints);
cosa = cosd(angle); sina = sind(angle);
rotMats = zeros(2, 2, npoints);
for i=1:npoints
    rotMats(:, :, i) = [ cosa(i)  sina(i);
                        -sina(i)  cosa(i)];
end
minwidth = inf; mini = NaN;
for i=1:npoints
    % Compute rotation matrix to orient the convex hull with the main x-y
    % axis (this way it is easier to find the bounding box with the min and
    % max values of x/y)
    rotVertices = (rotMats(:, :, i)*(vertices - [cx, cy])')';

    % Obtain width length of the rotated polygon
    maxx = max(rotVertices(:, 1)); minx = min(rotVertices(:, 1));
    l = maxx - minx;
    if l < minwidth
        minwidth = l;
        mini = i;
    end
end

% Return minimum width direction
thetamin = angle(mini);
height = area(polyshape(x, y))/minwidth;
% Constrain angle between 0º and 180º
if thetamin >= 180
    thetamin = thetamin - 180;
end
axes = [cosd(thetamin), sind(thetamin)];

end