function varargout = interppolygon(varargin)
% Given a set of polygon vertices, this function performs a linear
% interpolation of those to fulfill any possible gap in the polygon
% boundary (more density).
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        [x, y] = interppolygon(x, y)
%               [x, y, z] = interppolygon(x, y, z)
%
% Inputs:
%   > x:        If nargin > 2, this is the rectangular x-coordinates array
%               of the polygon boundary vertices. Otherwise, this 
%               corresponds to the longitude values array of the polygon 
%               boundary vertices 
%   > y:        If nargin > 2, this is the rectangular y-coordinates array
%               of the polygon boundary vertices. Otherwise, this 
%               corresponds to the latitude values array of the polygon 
%               boundary vertices 
%   > z:        rectangular z-coordinates of the polygon boundary vertices
%
% Outputs:
%   > x:        updated x-coordinates or longitude values array of the 
%               polygon boundary vertices (interpolated)
%   > y:        updated y-coordinates or latitude values array of the 
%               polygon boundary vertices (interpolated)
%   > z:        updated z-coordinates of the polygon boundary vertices
%               (interpolated)

% Variable allocation and error check
x = varargin{1};
y = varargin{2};
if length(x) ~= length(y)
    error("Input arguments must have the same length");
end
if nargin == 3
    z = varargin{3};
    if length(z) ~= length(x)
        error("Input arguments must have the same length");
    end
elseif nargin > 3
    error("Too many input arguments");
end

% Sort the boundary vertices in clockwise order (ccw would also work) to
% enable the linear interpolation
if nargin > 2
    [x, y, z] = sortcw(x, y, z);
    vertices = [x y z];
else
    [x, y] = sortcw(x, y);
    vertices = [x y];
end
vertices(end+1, :) = vertices(1, :); % close polygon

% Boundary mesh density (resolution)
dist = zeros(length(vertices), 1);
for i=2:length(vertices)
    dist(i) = norm(vertices(i, :) - vertices(i-1, :));
end
maxdiff = mean(dist);

% Interpolate polygon boundary vertices
aux = [];
for i=1:length(vertices)-1
    v = vertices(i+1, :) - vertices(i, :); % director vector
    if norm(v) > 2*maxdiff
        % Linear (approximation) interpolation between vertices
        lambda = linspace(0, 1, ceil(norm(v)/maxdiff)); % line parametrization
        for l=1:length(lambda)
            aux(end+1, :) = vertices(i, :) + v.*lambda(l);
        end
    end
end

% Add all the interpolated vertices to the coordinate matrix
vertices((end+1):(length(vertices)+length(aux)), :) = aux;

% Sort the polygon boundary vertices (this is mainly for representation)
if nargin > 2
    [varargout{1}, varargout{2}, varargout{3}] = sortcw(vertices(:,1), ...
        vertices(:,2), vertices(:,3));
else
    [varargout{1}, varargout{2}] = sortcw(vertices(:,1), vertices(:,2));
end

end