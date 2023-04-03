function [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma0, sg0, gamma, targetArea, gridPoints0, gridPoints, method)
% Flood-fill recursive algorithm that discretizes the target area by
% "flooding" the region with 2D rectangular elements. The grid is 
% determined by the input width, height and overlaps in both directions.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% 
% Usage:        gridPoints = floodFillVoronoiAlgorithm(w, h, ovlapx, ovlapy,
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
visited= false;

%
fp = footprint(gamma(1), gamma(2), t, inst, sc, target, theta);
w = sx0;
h = sy0;

ovlapx = olapx*w/100; ovlapy = olapy*h/100; % convert overlaps from 
% percentage to degrees of latitude and longitude, respectively

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
    gridPoints0(end+1, :) = gamma;
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
    % Check if the cell has been previously visited
    for i=1:size(gridPoints0,1)
        if norm(gridPoints0(i,:) - gamma) < 1e-5
            visited = true;
            break;
        end
    end

    % Disregard those cases where the footprint does not cover a certain
    % minimum of the target (this also avoids sub-optimality in the
    % optimization algorithms)
    areaInter = areaT - areaI;
    fpArea = area(fpshape);

    if visited || areaInter/fpArea < 0.1
        return;
    end

    % In case it has not been previously visited, then check the cardinal
    % (and diagonal neighbors in case the method is set to 8fill)
    % neighbors recursively
    s1 = sg0;
    if isequal(gamma0(1), gamma(1))
        gammap = gamma;
        s2 = sx0;
        s2y = sy0;
    else
        s2 = fp.sizex;
        s2y = fp.sizey;
        gammaold = gridPoints(end, :);
        if gamma(1) > gammaold(1)
            deltax = gamma(1) - gammaold(1) - .5*(s1 + s2);
            gammap(1) = gamma(1) - deltax;
        else
            deltax = -(gamma(1) - gammaold(1)) - .5*(s1 + s2);
            gammap(1) = gamma(1) + deltax;
        end
        gammap(2) = gamma(2);
    end
    gridPoints0(end+1,:) = gamma;
    gridPoints(end+1,:) = gammap;
      pp = polyshape([gammap(1)-s2/2, gammap(1)-s2/2, gammap(1)+s2/2, gammap(1) + s2/2, gammap(1)-s2/2],[gammap(2)+ s2y/2, gammap(2)- s2y/2, gammap(2)- s2y/2, gammap(2)+ s2y/2, gammap(2)+ s2y/2]);
      plot(pp, 'FaceColor', [0.93,0.69,0.13], 'FaceAlpha', 0.2);
      plot(gamma(1), gamma(2), 'b*')
      plot(gammap(1), gammap(2), 'r^')
      drawnow

    [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1)-w+ovlapx,          gamma(2)], targetArea, gridPoints0, gridPoints, method); % west
    [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1),          gamma(2)-h+ovlapy], targetArea, gridPoints0, gridPoints, method); % south
        [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1)+w-ovlapx,          gamma(2)], targetArea, gridPoints0, gridPoints, method); % east
    [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1),          gamma(2)+h-ovlapy], targetArea, gridPoints0, gridPoints, method); % north
    if isequal(method,'8fill')
        [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1)-w+ovlapx, gamma(2)+h-ovlapy], targetArea, gridPoints0, gridPoints, method); % northwest
        [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1)-w+ovlapx, gamma(2)-h+ovlapy], targetArea, gridPoints0, gridPoints, method); % southwest
        [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1)+w-ovlapx, gamma(2)+h-ovlapy], targetArea, gridPoints0, gridPoints, method); % northeast
        [gridPoints0, gridPoints] = floodFillVoronoiAlgorithm(sx0, sy0, t, inst, sc, target,...
    theta, olapx, olapy, gamma, s2, [gamma(1)+w-ovlapx, gamma(2)-h+ovlapy], targetArea, gridPoints0, gridPoints, method); % southeast
    end
end

% if isempty(gridPoints)
%     error("Empty grid!!")
% end
end