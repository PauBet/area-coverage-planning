function bbox = smallestBoundingBox(varargin)
% Given a set of polygon vertices, this function computes the smallest
% bounding box that encloses the set of bounding vertices
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        bbox = smallestBoundingBox(x, y)
%
% Inputs:
%   > x:        x coordinates of the polygon vertices (in our case,
%               longitude)
%   > y:        y coordinates of the polygon vertices (in our case,
%               latitude)
%   > angle:    in case the 
%
% Outputs:
%   > bbox:     struct containing main parameters of the smallest bounding
%               box
%       # size1: rectangle edge 1
%       # size2: rectangle edge 2
%       # boundPoints: boundary points that define the rectangle perimeter

%
x = varargin{1};
y = varargin{2};
if nargin > 2
    angle = deg2rad(varargin{3});
else
    angle = NaN;
end

% Check if the polygon is divided in two (a.m. intersection)...
ind = find(isnan(x));
if ~isempty(ind)
    x(1:ind) = x(1:ind) + 360;
    x(ind) = [];
    y(ind) = [];
end

% Compute the convex hull of the polygon vertices (fortunately, matlab has
% a function that computes it)
k = convhull(x, y);

% Save the vertices of the convex hull and find its centroid
vertices = [x(k), y(k)];
[cx, cy] = centroid(polyshape(vertices(:,1), vertices(:,2)));

% Sort the vertices in clockwise direction
[vertices(:,1), vertices(:,2)] = sortcw(vertices(:,1), vertices(:,2));
vertices(end+1, :) = vertices(1, :); % close polygon

% Smallest surrounding rectangle
rotVertices = zeros(length(vertices), 2);
minArea = inf;
if isnan(angle)
    for i=2:length(vertices)
        % For each edge of the convex hull
        edge = vertices(i, :) - vertices(i-1, :);

        % Compute edge orientation
        angle = atan2(edge(2), edge(1));

        % Compute rotation matrix to orient the convex hull with the main x-y
        % axis (this way it is easier to find the bounding box with the min and
        % max values of x/y)
        rotmat = [cos(angle) sin(angle);
            -sin(angle)  cos(angle)];
        for j=1:length(vertices)
            rotVertices(j, :) = rotmat*(vertices(j, :)' - [cx, cy]');
        end

        % Obtain the bounding box
        maxx = max(rotVertices(:, 1)); minx = min(rotVertices(:, 1));
        maxy = max(rotVertices(:, 2)); miny = min(rotVertices(:, 2));
        boundingArea = (maxx - minx)*(maxy - miny);

        % Compare and find the minimum bounding box area, and save the
        % parameters in the output struct bbox
        if boundingArea < minArea
            minArea = boundingArea;
            boundPoints = [maxx maxy;
                maxx miny;
                minx miny;
                minx maxy];
            minAngle = angle;
            opt_maxx = maxx;
            opt_maxy = maxy;
            opt_minx = minx;
            opt_miny = miny;
        end
    end
else
    % Compute rotation matrix to orient the convex hull with the main x-y
    % axis (this way it is easier to find the bounding box with the min and
    % max values of x/y)
    rotmat = [cos(angle) sin(angle);
        -sin(angle)  cos(angle)];
    for j=1:length(vertices)
        rotVertices(j, :) = rotmat*(vertices(j, :)' - [cx, cy]');
    end

    % Obtain the bounding box
    maxx = max(rotVertices(:, 1)); minx = min(rotVertices(:, 1));
    maxy = max(rotVertices(:, 2)); miny = min(rotVertices(:, 2));
    boundPoints = [maxx maxy;
        maxx miny;
        minx miny;
        minx maxy];
    opt_maxx = maxx;
    opt_maxy = maxy;
    opt_minx = minx;
    opt_miny = miny;
    minAngle = angle;
end

% Save bound points of the smallest bounding box
rotmat = [cos(minAngle) sin(minAngle);
         -sin(minAngle) cos(minAngle)];
rotBoundPoints = zeros(length(boundPoints) + 1, 2);
for bb=1:length(boundPoints)
    rotBoundPoints(bb,:) = transpose(rotmat)*boundPoints(bb, :)'...
        + [cx, cy]';

    % Check a.m. intercept case
    if ~isempty(ind) && rotBoundPoints(bb, 1) >= 180
        rotBoundPoints(bb, 1) = rotBoundPoints(bb, 1) - 360;
    end
end
rotBoundPoints(end, :) = rotBoundPoints(1, :); % close polygon
[bbox.boundPoints(:,1), bbox.boundPoints(:,2)] = ...
    sortcw(rotBoundPoints(:,1), rotBoundPoints(:,2)); % sort
% the vertices in clockwise direction

% Find the footprint's x-axis angle w.r.t. longitude axis
if nargin < 3
    edge = [opt_maxx 0];
    % if (edge(1) < 0 && edge(2) < 0) || (edge(1) > 0 && edge(2) < 0)
    %     edge = -edge;
    % end
    angle = acos(dot([1 0], edge) / norm(edge));
end

% Save min/max longitude and latitude values
bbox.maxlon = opt_maxx;
bbox.minlon = opt_minx;
bbox.maxlat = opt_maxy;
bbox.minlat = opt_miny;

% Save smallest bounding box size
bbox.size1 = opt_maxx - opt_minx; bbox.size2 = opt_maxy - opt_miny;

% Save bounding box orientation
bbox.angle = rad2deg(angle);

end