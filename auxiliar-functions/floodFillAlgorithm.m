function gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, gamma,...
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
%                   Units are irrelevant as long as they are consistent
%   > ovlapy:       grid gootprint overlap in the vertical direction. Units
%                   are irrelevant as long as they are consistent
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
visited= false;

% Rectangular element definition
fpx = [gamma(1) - w/2, gamma(1) - w/2, gamma(1) + w/2, gamma(1) + w/2];
fpy = [gamma(2) + h/2, gamma(2) - h/2, gamma(2) - h/2, gamma(2) + h/2];

% Subtract the rectangle from the targetArea
inter = subtract(polyshape(targetArea(:,1), targetArea(:,2)), ...
    polyshape(fpx, fpy));
areaI = area(inter);
targetpshape = polyshape(targetArea(:,1), targetArea(:,2));
areaT = area(targetpshape);

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
    % Check if the cell has been previously visited
    for i=1:size(gridPoints,1)
        if abs(gridPoints(i,:) - gamma) < 1e-5
            visited = true;
            break;
        end
    end

    % Disregard those cases where the footprint does not cover a certain
    % minimum of the target (this also avoids sub-optimality in the
    % optimization algorithms)
    areaInter = areaT - areaI;
    fpArea = area(polyshape(fpx, fpy));

    % In case it has not been previously visited, then check the cardinal
    % (and diagonal neighbors in case the method is set to 8fill) 
    % neighbors recursively
    if ~visited && areaInter/fpArea >= 0.2
        gridPoints(end+1,:) = gamma;
        plot([gamma(1)-w/2, gamma(1)-w/2, gamma(1)+w/2, gamma(1) + w/2, gamma(1)-w/2],[gamma(2)+ h/2, gamma(2)- h/2, gamma(2)- h/2, gamma(2)+ h/2, gamma(2)+ h/2],'Color','g');
        drawnow

        gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1)-w+ovlapx,          gamma(2)], targetArea, gridPoints, method); % west
        gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1),          gamma(2)-h+ovlapy], targetArea, gridPoints, method); % south
        gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1),          gamma(2)+h-ovlapy], targetArea, gridPoints, method); % north
        gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1)+w-ovlapx,          gamma(2)], targetArea, gridPoints, method); % east
        if isequal(method,'8fill')
            gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1)-w+ovlapx, gamma(2)+h-ovlapy], targetArea, gridPoints, method); % northwest
            gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1)-w+ovlapx, gamma(2)-h+ovlapy], targetArea, gridPoints, method); % southwest
            gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1)+w-ovlapx, gamma(2)+h-ovlapy], targetArea, gridPoints, method); % northeast
            gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, [gamma(1)+w-ovlapx, gamma(2)-h+ovlapy], targetArea, gridPoints, method); % southeast
        end
    end
end

if isempty(gridPoints)
    error("Empty grid!!")
end
end