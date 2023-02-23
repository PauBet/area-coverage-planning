function bbox = smallestBoundingBox(x, y)
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
%
% Outputs:
%   > bbox:     struct containing main parameters of the smallest bounding
%               box
%       # size1: rectangle edge 1
%       # size2: rectangle edge 2
%       # boundPoints: boundary points that define the rectangle perimeter

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

% Save the vertices of the convex hull
vertices(:, 1) = x(k); vertices(:, 2) = y(k);

% Sort the vertices in clockwise direction
[vertices(:,1), vertices(:,2)] = sortcw(vertices(:,1), vertices(:,2));
vertices(end+1, :) = vertices(1, :); % close polygon

% Smallest surrounding rectangle
rotVertices = zeros(length(vertices), 2);
minArea = inf;
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
    [cx, cy] = centroid(polyshape(vertices(:,1), vertices(:,2)));
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

        % Save bound points of the smallest bounding box
        for bb=1:length(boundPoints)
            bbox.boundPoints(bb,:) = transpose(rotmat)*boundPoints(bb, :)'...
                + [cx, cy]';

            % Check a.m. intercept case
            if ~isempty(ind) && bbox.boundPoints(bb, 1) >= 180 
                bbox.boundPoints(bb, 1) = bbox.boundPoints(bb, 1) - 360;
            end
        end
        bbox.boundPoints(end+1, :) = bbox.boundPoints(1,:); % close polygon
        [bbox.boundPoints(:,1), bbox.boundPoints(:,2)] = ...
            sortcw(bbox.boundPoints(:,1), bbox.boundPoints(:,2)); % sort
            % the vertices in clockwise direction
        
        % Save min/max longitude and latitude values
        bbox.maxlon = maxx;
        bbox.minlon = minx;
        bbox.maxlat = maxy;
        bbox.minlat = miny;

        % Save smallest bounding box size
        bbox.size1 = maxx - minx; bbox.size2 = maxy - miny;

        % Save angle
        if angle >= pi
            angle = angle - pi;
        elseif angle < 0
            angle = angle + pi;
        end
        
        % first quadrant
        if angle > pi/2 && angle <= pi
            angle = angle - pi/2;
        end
        bbox.angle = rad2deg(angle);
    end
end

end