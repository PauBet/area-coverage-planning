function [gridPoints, dirx, diry] = grid2Dbis(fpref, ovlapx, ovlapy, ...
    gamma, targetArea)
% Grid discretization (using flood-fill algorithm) of a region of interest
% given a reference footprint (unit measure to create the allocatable
% cells)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% 
% Usage:        matrixGrid = grid2D(fpref, ovlapx, ovlapy, gamma,...
%                   targetArea)
%
% Inputs:
%   > fpref:        struct containing the parameters that define the 
%                   footprint. In this function, the following are needed:
%      # sizex:     footprint size in the x direction (longitude), in deg
%      # sizey:     footprint size in the y direction (latitude), in deg
%                   See function 'footprint' for further information.
%   > olapx:        grid footprint overlap in the x direction, in 
%                   percentage
%   > olapy:        grid footprint overlap in the y direction, in 
%                   percentage
%   > targetArea:   matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D.
%       # targetArea(:,1) correspond to the x values of the vertices
%       # targetArea(:,2) correspond to the y values of the vertices
% 
% Outputs:
%   > matrixGrid:   cell matrix containing the grid discretization of the
%                   region-of-interest (ROI).
%                   Each point is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
% The matrix sorts the discretized points (flood-fill) by latitude and
% longitude according to the following structure:
%
%               longitude
%               (-) --------> (+)                   
%  latitude (+) [a11]  [a12] ⋯
%            ¦  [a21]
%            ¦    ⋮
%            ∨
%           (-)   

% Pre-allocate variables
gridPoints = {};

% Get the footprint angle, i.e., the angle that the 2D footprint forms with
% respect to the meridian-equator axes
angle = deg2rad(-fpref.angle);

% Filling the region-of-interest (roi) with a footprint that is not aligned
% with the meridian-equator axes is equivalent to filling the oriented
% target area with an aligned footprint (angle = 0). Therefore, we rotate
% the region-of-interest to orient it according to the footprint
rotmat = [cos(angle)   -sin(angle);
          sin(angle)   cos(angle)];
% matrixGrid directions x and y
dirx = rotmat(1, :);
diry = rotmat(2, :);
[cx, cy] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
orientedArea  = zeros(length(targetArea), 2);
for j=1:length(targetArea)
    orientedArea(j, :) = [cx, cy]' + rotmat*(targetArea(j, :)' - ...
        [cx, cy]');
end
gamma = [cx, cy]' + rotmat*(gamma' - [cx, cy]');

% If the area is divided in smaller regions, then we get the convex polygon
% that encloses all of them (flood-fill)
aux(:, 1) = orientedArea(~isnan(orientedArea(:, 1)), 1); % convhull 
% does not accept NaN nor Inf
aux(:, 2) = orientedArea(~isnan(orientedArea(:, 2)), 2);
k = convhull(aux(:, 1), aux(:, 2)); % boundary vertices that constitute the
% convex polygon
periArea(:, 1) = aux(k, 1);
periArea(:, 2) = aux(k, 2);

% Auxiliary figure
% figure
% hold on;
% plot(polyshape(orientedArea(:,1), orientedArea(:,2)))
% plot(gamma(1), gamma(2), 'r*')
% plot(polyshape(periArea(:, 1), periArea(:, 2)), 'FaceColor', 'none')
% drawnow

% Flood-fill algorithm to get the grid points of the oriented roi
% gridPoints = floodFillAlgorithmPar(fpref.sizex, fpref.sizey, ovlapx, ...
%  ovlapy, gamma, orientedArea, gridPoints, '8fill');
gridPoints = floodFillAlgorithm(fpref.width, fpref.height, ovlapx, ...
 ovlapy, gamma, orientedArea, periArea, [], [], '4fill');

% Convert into topographical coordinates
grid = inst2topo(gridPoints, rotmat, target, et, sc);

% Sort matrix


end