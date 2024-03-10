function gridPoints = floodFillAlgorithmPar(w, h, olapx, olapy, gamma,...
    targetArea, gridPoints, method)
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
%   > ovlapx:       grid footprint overlap in the horizontal direction.
%                   Units are in percentage of width
%   > ovlapy:       grid gootprint overlap in the vertical direction. Units
%                   are in percentage of height
%   > gamma:        grid origin point
%   > targetArea:   matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D 
%       # targetArea(:,1) correspond to the x values of the vertices
%       # targetArea(:,2) correspond to the y values of the vertices
%   > gridPoints:   matrix containing the discretized grid points of the
%                   region-of-interest. The recursive calls of the
%                   algorithm will fill this matrix. These grid points
%                   represent the center of the rectangular elements used
%                   to fill the region
%   > method:       string name of the method. '4fill' fills the roi by
%                   searching the cardinal directions. '8fill' consider
%                   also the diagonal neighbours
% 
% Outputs:
%   > gridPoints:   matrix containing the discretized gridPoints of the
%                   region-of-interest. The recursive calls of the
%                   algorithm will fill this matrix. These grid points
%                   represent the center of the rectangular elements used
%                   to fill the region

% Variables pre-allocation
inside = false;
ovlapx = olapx*w/100; ovlapy = olapy*h/100; % convert overlaps from 
% percentage to degrees of latitude and longitude, respectively

% Check if the cell has been previously visited
for i=1:size(gridPoints,1)
    if abs(gridPoints(i,:) - gamma) < 1e-5
        return;
    end
end

% Rectangular element definition
fpx = [gamma(1) - w/2, gamma(1) - w/2, gamma(1) + w/2, gamma(1) + w/2];
fpy = [gamma(2) + h/2, gamma(2) - h/2, gamma(2) - h/2, gamma(2) + h/2];

% Subtract the rectangle from the targetArea
targetpshape = polyshape(targetArea(:,1), targetArea(:,2));
fpshape = polyshape(fpx, fpy);
inter = subtract(targetpshape, fpshape);
areaI = area(inter);
areaT = area(targetpshape);

% Check: the footprint is larger than the region of interest...
if areaI == 0
    gridPoints(end+1, :) = gamma;
    return;
end

% Check if the rectangle at gamma and size [w,h] is contained in
% the target area (either partially or totally)
if inpolygon(gamma(1), gamma(2), targetArea(:,1), targetArea(:,2))
    inside = true;
else
    if abs(areaI - areaT) > 1e-5
            inside = true;
    end
end

if inside
    % Disregard those cases where the footprint does not cover a certain
    % minimum of the target (this also avoids sub-optimality in the
    % optimization algorithms)
    areaInter = areaT - areaI;
    fpArea = area(fpshape);

    if areaInter/fpArea < 0.2
        return;
    end

    % In case it has not been previously visited, then check the cardinal
    % (and diagonal neighbors in case the method is set to 8fill)
    % neighbors recursively
    gridPoints(end+1, :) = gamma;

    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1)-w+ovlapx,          gamma(2)], targetArea, gridPoints, method); % west
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1),          gamma(2)-h+ovlapy], targetArea, gridPoints, method); % south
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1),          gamma(2)+h-ovlapy], targetArea, gridPoints, method); % north
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1)+w-ovlapx,          gamma(2)], targetArea, gridPoints, method); % east
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1)-w+ovlapx, gamma(2)+h-ovlapy], targetArea, gridPoints, method); % northwest
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1)-w+ovlapx, gamma(2)-h+ovlapy], targetArea, gridPoints, method); % southwest
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1)+w-ovlapx, gamma(2)+h-ovlapy], targetArea, gridPoints, method); % northeast
    gridPoints = fetchOutputs(F);
    F = parfeval(@floodFillAlgorithmPar, 1, w, h, olapx, olapy, [gamma(1)+w-ovlapx, gamma(2)-h+ovlapy], targetArea, gridPoints, method); % southeast
    gridPoints = fetchOutputs(F);
end
end