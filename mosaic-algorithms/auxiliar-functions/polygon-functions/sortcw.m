function varargout = sortcw(varargin)
% Given a set of polygon vertices, this function sorts them in clockwise
% order
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        [x, y] = sortcw(x, y)
%               [x, y, z] = sortcw(x, y, z)
%
% Inputs:
%   > x:        x rectangular coordinate (nargin = 3) or longitude values
%               (nargin = 2). Units are irrelevant as long as the 2 or 3
%               arrays are consistent
%   > y:        y rectangular coordinate (nargin = 3) or latitude values
%               (nargin = 2). Units are irrelevant as long as the 2 or 3
%               arrays are consistent
%   > z:        z rectangular coordinate (nargin = 3)
%
% Outputs:
%   > Same as inputs but sorted in clockwise order
%
% Note: 2D sorting algorithm does not work with concave algorithms. In such
% case, 3D rectangular coordinates are recommended.
%
if nargin < 3
    % This algorithm consists of dividing the space in 4 quadrants,
    % centered at the polygon centroid (any other inner point would be also
    % valid). We calculate iteratively the angle with respect to that point
    % and sort the vertices according to their respective angle value.
    % This algorithm does not work with non-convex polygons
    x = varargin{1}; y = varargin{2};
    [cx, cy] = centroid(polyshape(x, y)); % polygon centroid
    angle = atan2(y - cy, x - cx); % obtain angle
    [~, ind] = sort(angle, 'descend'); % sort angles and get the indices
    x = x(ind); % sort the longitude values according to the angle order
    y = y(ind); % sort the latitude values according to the angle order
    varargout{1} = x; varargout{2} = y; % save output
    
elseif nargin == 3
    % Algorithm extracted from 
    % https://stackoverflow.com/questions/47949485/
    % sorting-a-list-of-3d-points-in-clockwise-order
    x = varargin{1}; y = varargin{2}; z = varargin{3};
    innerpoint = [mean(x), mean(y), mean(z)]; % find an inner point (this 
    % is not the centroid but still works)
    i = [1, 0, 0]; j = [0, 1, 0]; k = [0, 0, 1];

    % Compute the cross products (perpendicular to the surface normal)
    pn(:, 1) = cross(i, innerpoint); pn(:, 2) = cross(j, innerpoint); 
    pn(:, 3) = cross(k, innerpoint);
    
    % Take the largest cross product
    [~, pind] = max([norm(pn(:, 1)), norm(pn(:, 2)), norm(pn(:, 3))]);
    p = pn(:, pind);

    % Compute vector perpendicular both the surface normal and p
    q = cross(innerpoint, p);

    % Now we have two perpendicular reference vectors in the plane given by
    % the normal. Take triple products of those, and these will be the sine
    % and the cosine of an angle that can be used for sorting
    angle = zeros(1, length(x));
    for ii=1:length(x)
        rmc = [x(ii), y(ii), z(ii)] - innerpoint;
        t = dot(innerpoint,(cross(rmc, p)));
        u = dot(innerpoint,(cross(rmc, q)));
        angle(ii) = atan2(u, t);
    end

    % Sort angles and get the indices to sort the vertices
    [~, ind] = sort(angle);
    x = x(ind); y = y(ind); z = z(ind);
    varargout{1} = x; varargout{2} = y; varargout{3} = z;
else
    error("Too many input arguments")
end
end