function gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, gamma, targetArea, gridPoints, method)
%% Variable initialization
inside = false;
visited= false;
% Footprint definition
fpx = [gamma(1) - w/2, gamma(1) - w/2, gamma(1) + w/2, gamma(1) + w/2];
fpy = [gamma(2) + h/2, gamma(2) - h/2, gamma(2) - h/2, gamma(2) + h/2];
% Subtract the footprint from the targetArea
inter = subtract(polyshape(targetArea(:,1), targetArea(:,2)), polyshape(fpx, fpy));
areaI = area(inter);
targetpshape = polyshape(targetArea(:,1), targetArea(:,2));
areaT = area(targetpshape);

%% Check if the footprint at gamma and size [w,h] is contained in
%% the target area (either partially or totally)
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
    % (and diagonal neighbors in case the method is set to 8fill) recursively
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