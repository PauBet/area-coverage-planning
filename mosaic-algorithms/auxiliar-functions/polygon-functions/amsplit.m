function [xf, yf] = amsplit(x, y)
% Provided the vertices of a polygon in latitudinal coordinates, this 
% function analyzes if the polygon intercepts the anti-meridian line and, 
% in that case, divides the polygon by this line
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% 
% Usage:        [xf, yf] = amsplit(x, y)
%
% Inputs:
%   > x:        array of longitude values in [deg]. x ∈ [-180, 180]
%   > y:        array of latitude values in [deg]. y ∈ [-90, 90]
%
% Outputs:
%   > x:        array of longitude values in [deg]. In case the polygon 
%               intercepts the a.m. line, then the longitude values of the 
%               polygon are separated by a NaN value. x ∈ [-180, 180]
%   > y:        array of latitude values in [deg]. In case the polygon 
%               intercepts the a.m. line, then the latitude values of the 
%               polygon are separated by a NaN value. y ∈ [-90, 90]

% [To be resolved]: If polygon longitude size is 180º, the function does
% nothing.
if (max(x) - min(x)) == 180
    warning("Full longitude. This function does nothing, in this case." + ...
        "The user must check if the polygon is well defined!");
end

% Previous check...
if (max(x) - min(x)) <= 180
    xf = x;
    yf = y;
    return;
end

% In case the polygon indeed intercepts this line, then we're going to
% calculate the intercept points by computing the intersection between two
% polygons: the input polygon and another one which corresponds to the
% anti-meridian line (actually, it's a small polygon because the MATLAB's
% 'intersect' function does not operate well with lines)
x(x < 0) = x(x < 0) + 360; % continuous polygon
poly1 = polyshape(x, y);
vpoly2(:, 1) = [180*ones(1, 20) 181*ones(1, 20)];
vpoly2(:, 2) = [linspace(-90, 90, 20) linspace(90, -90, 20)];
poly2 = polyshape(vpoly2(:, 1), vpoly2(:, 2));

% Compute the intersection points
polyinter = intersect(poly1, poly2);
xinter = polyinter.Vertices(:, 1);
yinter = polyinter.Vertices(:, 2);

% Only keep the anti-meridian intercepts
yi = yinter(abs(xinter - 180) < 1e-2);
xi = 180*ones(length(yi), 1);

% Define the new polygon
P = [x y];

% Add the intersection points to the polygon vertices
P = [P; xi yi];

% Split the polygon in two, cleaved by the anti-meridian line
x = sort(P(:, 1));
P = sortrows(P);
poly1 = P(x >= 180, :); % this is the polygon that falls in negative lon.
poly2 = P(x <= 180, :); % this is the polygon that falls in positive lon.

% Sort the polygon vertices in clockwise order...
[poly1(:,1), poly1(:,2)] = sortcw(poly1(:,1), poly1(:,2));
[poly2(:,1), poly2(:,2)] = sortcw(poly2(:,1), poly2(:,2));

% Retrieve the original values (longitude values cannot be >180º in our
% system) and output the final vertices
xf = poly1(:, 1) - 360;
xf(end + 1) = NaN;
xf(end + 1:(length(xf) + length(poly2(:, 1)))) = poly2(:, 1);

yf = poly1(:, 2);
yf(end + 1) = NaN;
yf(end + 1:(length(yf) + length(poly2(:, 2)))) = poly2(:, 2);
end