function [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, gamma,...
    targetArea, perimeterArea, gridPoints, vPoints, method)
% Flood-fill recursive algorithm that discretizes the target area by
% "flooding" the region with 2D rectangular elements. The grid is 
% determined by the input width, height and overlaps in both directions.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% 
% Usage:        gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy,
%                   gamma, targetArea, gridPoints, method)
%
% Inputs:
%   > w:            horizontal resolution. Units are irrelevant as long as
%                   they are consistent
%   > h:            vertical resolution. Units are irrelevant as long as
%                   they are consistent
%   > olapx:        grid footprint overlap in the horizontal direction.
%                   Units are in percentage of width
%   > olapy:        grid gootprint overlap in the vertical direction. Units
%                   are in percentage of height
%   > gamma:        grid origin point (seed)
%   > targetArea:   matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D 
%       # targetArea(:,1) correspond to the x values of the vertices
%       # targetArea(:,2) correspond to the y values of the vertices
%   > perimeterArea:matrix containing the vertices of the polygon that
%                   encloses all of the uncovered area. At the beginning:
%                   perimeterArea = targetArea, but as the observations
%                   advance, this is going to change. Recommended: use the
%                   convex hull function of the uncovered area.
%   > gridPoints:   matrix containing the discretized grid points of the
%                   region-of-interest. The recursive calls of the
%                   algorithm will fill this matrix. These grid points
%                   represent the center of the rectangular elements used
%                   to fill the region.
%        # When calling this function: gridPoints = [];
%   > vPoints:      matrix containing the visited points (to prevent
%                   gridlock)
%        # When calling this function: vPoints = [];
%   > method:       string name of the method. '4fill' fills the roi by
%                   searching the cardinal directions. '8fill' considers
%                   also the diagonal neighbours
% 
% Outputs:
%   > gridPoints:   matrix containing the discretized gridPoints of the
%                   region-of-interest. The recursive calls of the
%                   algorithm will fill this matrix. These grid points
%                   represent the center of the rectangular elements used
%                   to fill the region
%   > vPoints:      matrix containing the visited points (to prevent
%                   gridlock)
%
% Note: this function creates a convex polygon that encloses all the
% uncovered area (even when this is divided in portions) and tours the
% whole area. It is less computationally efficient than the classic
% flood-fill algorithm, but it is convenient to prevent sub-optimal
% fillings of the uncovered area (isolated points).

% Variables pre-allocation
inside = false;
ovlapx = olapx*w/100; ovlapy = olapy*h/100; % convert overlaps from 
% percentage to degrees of latitude and longitude, respectively
epsilon = 0.05;

% Check if the cell has been previously visited
for i=1:size(vPoints, 1)
    if norm(vPoints(i,:) - gamma) < 1e-5
        return;
    end
end
% otherwise...
vPoints(end+1, :) = gamma;

% Rectangular element definition
fpx = [gamma(1) - w/2, gamma(1) - w/2, gamma(1) + w/2, gamma(1) + w/2];
fpy = [gamma(2) + h/2, gamma(2) - h/2, gamma(2) - h/2, gamma(2) + h/2];

% Subtract the allocated cell (footprint) from the perimeterArea
peripshape = polyshape(perimeterArea(:, 1), perimeterArea(:, 2));
fpshape = polyshape(fpx, fpy);
inter = subtract(peripshape, fpshape);
areaI = area(inter);
areaP = area(peripshape);

% Check: the footprint is larger than the region of interest...
if areaI == 0
    gridPoints(end+1, :) = gamma;
    return;
end

% Check if the rectangle at gamma and size [w,h] is contained in
% the perimeter area (either partially or totally)
if inpolygon(gamma(1), gamma(2), targetArea(:,1), targetArea(:,2))
    inside = true;
else
    if abs(areaI - areaP)/area(fpshape) > epsilon
            inside = true;
    end
end

if inside
    % Disregard those cases where the footprint does not cover a certain
    % minimum of the roi (this also avoids sub-optimality in the
    % optimization algorithms)
    targetpshape = polyshape(targetArea(:, 1), targetArea(:, 2));
    areaT = area(targetpshape);
    inter = subtract(targetpshape, fpshape);
    areaI = area(inter);
    areaInter = areaT - areaI;
    fpArea = area(fpshape);
    
    if areaInter/fpArea > epsilon
        gridPoints(end+1, :) = gamma;
        % pp = polyshape([gamma(1)-w/2, gamma(1)-w/2, gamma(1)+w/2, gamma(1) + w/2, gamma(1)-w/2],[gamma(2)+ h/2, gamma(2)- h/2, gamma(2)- h/2, gamma(2)+ h/2, gamma(2)+ h/2]);
        % plot(pp, 'FaceColor', [0.93,0.69,0.13], 'FaceAlpha', 0.2);
        % plot(gamma(1), gamma(2), 'r^')
        % drawnow
    else
        if isempty(gridPoints)
            return; % first iteration: the footprint located at the 
            % original gamma does not cover the minimum required
        end
        % plot(gamma(1), gamma(2), 'b^')
        % drawnow
    end

    % Check the cardinal (and diagonal neighbors in case the method is set
    % to 8fill) neighbors recursively
    
    [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1)-w+ovlapx,          gamma(2)], targetArea, perimeterArea, gridPoints, vPoints, method); % west
    [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1),          gamma(2)-h+ovlapy], targetArea, perimeterArea, gridPoints, vPoints, method); % south
    [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1),          gamma(2)+h-ovlapy], targetArea, perimeterArea, gridPoints, vPoints, method); % north
    [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1)+w-ovlapx,          gamma(2)], targetArea, perimeterArea, gridPoints, vPoints, method); % east
    if isequal(method,'8fill')
        [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1)-w+ovlapx, gamma(2)+h-ovlapy], targetArea, perimeterArea, gridPoints, vPoints, method); % northwest
        [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1)-w+ovlapx, gamma(2)-h+ovlapy], targetArea, perimeterArea, gridPoints, vPoints, method); % southwest
        [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1)+w-ovlapx, gamma(2)+h-ovlapy], targetArea, perimeterArea, gridPoints, vPoints, method); % northeast
        [gridPoints, vPoints] = floodFillAlgorithm(w, h, olapx, olapy, [gamma(1)+w-ovlapx, gamma(2)-h+ovlapy], targetArea, perimeterArea, gridPoints, vPoints, method); % southeast
    end
end
end