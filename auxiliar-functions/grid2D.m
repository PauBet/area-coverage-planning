function matrixGrid = grid2D(fpref, ovlapx, ovlapy, gamma, targetArea)
% Flood-fill grid discretization of the region-of-interest given a 
% footprint reference's size and orientation
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
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage
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

% Get the footprint bounding box
bbox  = smallestBoundingBox(fpref.bvertices(:, 1), fpref.bvertices(:, 2));

% Get the footprint angle, i.e., the angle that the 2D footprint forms with
% respect to the meridian-equator axes
angle = deg2rad(bbox.angle);

% Filling the region-of-interest (roi) with a footprint that is not aligned
% with the meridian-equator axes is equivalent to filling the oriented
% target area with an aligned footprint (angle = 0). Therefore, we rotate
% the region-of-interest to orient it according to the footprint
rotmat = [cos(angle)   -sin(angle);
          sin(angle)   cos(angle)];
[cx, cy] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
orientedArea  = zeros(length(targetArea), 2);
for j=1:length(targetArea)
    orientedArea(j, :) = [cx, cy]' + rotmat*(targetArea(j, :)' - ...
        [cx, cy]');
end
gamma = [cx, cy]' + rotmat*(gamma' - [cx, cy]');

% Flood-fill algorithm to get the grid points of the oriented roi
gridPoints = [];
gridPoints = floodFillAlgorithm(fpref.sizex, fpref.sizey, ovlapx, ...
    ovlapy, gamma, orientedArea, gridPoints, '4fill');

% Temp figure
% figure
% plot(polyshape(targetArea(:,1), targetArea(:,2)))
% hold on;
% plot(polyshape(orientedArea(:,1), orientedArea(:,2)))
% plot(gridPoints(:,1), gridPoints(:,2), 'b*')
% orientedGridPoints = zeros(length(gridPoints), 2);
% for j=1:length(gridPoints)
%     orientedGridPoints(j, :) = [cx, cy]' +  transpose(rotmat)*...
%         (gridPoints(j, :)' - [cx, cy]');
% end
% plot(orientedGridPoints(:,1), orientedGridPoints(:,2), 'r*')

% Build a matrix that will sort the gridPoints elements by latitude and
% longitude according to the following structure:
%
%               longitude
%               (-) --------> (+)                   
%  latitude (+) [a11]  [a12] ⋯
%            ¦  [a21]
%            ¦    ⋮
%            ∨
%           (-)             
sortedGrid = sortrows(gridPoints, -2); % the elements of gridPoints are
% sorted by latitude (+ to -)
uniqueLat = unique(sortedGrid(:,2)); % get the different latitude values
ind = abs(diff(uniqueLat)) < 1e-5; % double check that there are
% no "similar" latitude values (it may happen)
uniqueLat(ind) = [];
uniqueLon = unique(sortedGrid(:,1)); % get the different longitude values
% unique check
ind = abs(diff(uniqueLon)) < 1e-5; % double check that there are
% no "similar" longitude values (it may happen)
uniqueLon(ind) = [];

% Sort and rotate the grid points and insert them in the grid matrix
matrixGrid = cell(length(uniqueLat), length(uniqueLon));
for i=1:length(uniqueLat)
    % We will sweep across the grid by, first, latitude and, second,
    % longitude
    lat = uniqueLat(length(uniqueLat) + 1 - i);
    indlat = (abs(sortedGrid(:,2) - lat) < 1e-5);
    mrow = sort(sortedGrid(indlat));
    for j=1:length(mrow)
        indlon = abs(uniqueLon - mrow(j)) < 1e-5;
        lon = mrow(j);
        matrixGrid{i, indlon} = [cx, cy]' + transpose(rotmat)*...
            ([lon lat]' - [cx, cy]'); % rotate values to the original roi
            % orientation
    end
end

end