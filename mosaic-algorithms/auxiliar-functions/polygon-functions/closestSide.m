function [dir1, dir2] = closestSide(target, sc, t, targetArea, angle)
% Given a region-of-interest, this function defines what is the spacecraft
% ground track position with respect to the edges of the target area
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% 
% Usage:        cside = closestSide(target, sc, t, roi)
%
% Inputs:
%   > target:       SPICE string name of the target body
%   > sc:           SPICE string name of the spacecraft
%   > t:            time in TDB seconds past J2000 epoch
%   > roi:          matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D latitudinal coord. 
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices
% 
% Outputs:
%   > cside:        string defining the spacecraft ground track position
%                   with respect to the roi, i.e. 'up', 'down', 'left' or
%                   'right'

% Parameters
[~, targetFrame, ~] = cspice_cnmfrm(target); % target body-fixed reference 
% frame

% Rotate roi according to the footprint's angle
angle = -angle*cspice_rpd;
rotmat = [cos(angle)   -sin(angle);
          sin(angle)   cos(angle)];
[cx, cy] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
roi  = zeros(length(targetArea), 2);
for j=1:length(targetArea)
    roi(j, :) = [cx, cy]' + rotmat*(targetArea(j, :)' - ...
        [cx, cy]');
end

% Compute spacecraft sub-observer point to see what side of the target
% area is closer to it
% Assumption: Tri-axial ellipsoid to model the target surface
subobs = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t,...
    targetFrame, 'NONE', sc);
[~, sclon, sclat] = cspice_reclat(subobs); % latitudinal coordinates
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr; % [rad] to [deg]

aux = [cx, cy]' + rotmat*([sclon; sclat] - ...
        [cx, cy]');
sclon = aux(1); sclat = aux(2);

% Compute spacecraft sub-observer point afterwards to determine if the
% spacecraft is moving away or towards the closest side
% Assumption: Tri-axial ellipsoid to model the target surface
subobs_ = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t + 5*60,...
    targetFrame, 'NONE', sc);
[~, sclon_, sclat_] = cspice_reclat(subobs_); % latitudinal coordinates
sclon_ = sclon_*cspice_dpr; sclat_ = sclat_*cspice_dpr; % [rad] to [deg]

aux = [cx, cy]' + rotmat*([sclon_; sclat_] - ...
        [cx, cy]');
sclon_ = aux(1); sclat_ = aux(2);

% Find the 4 boundary vertices
maxlon = max(roi(:, 1)); minlon = min(roi(:, 1));
maxlat = max(roi(:, 2)); minlat = min(roi(:, 2));

% ROI's boundary box
xlimit = [minlon maxlon];
ylimit = [minlat  maxlat];
xbox = xlimit([1 1 2 2 1]);
ybox = ylimit([1 2 2 1 1]);

% Define line between the centroid and the ground track
x = [sclon cx];
y = [sclat cy];
[xi, yi] = polyxpoly(x, y, xbox, ybox);

if isempty(xi) % the ground track is inside the ROI's boundary box (no
    % intersection)

    p1 = [maxlon, maxlat];
    p2 = [maxlon, minlat];
    p3 = [minlon, minlat];
    p4 = [minlon, maxlat];

    % Calculate distance to the mid-points of the 4 edges of the boundary 
    % box
    midp(1, :) = .5*(p1 + p4); % midup
    midp(2, :) = .5*(p2 + p3); % middown
    midp(3, :) = .5*(p3 + p4); % midleft
    midp(4, :) = .5*(p1 + p2); % midright

    % Find minimum distance
    mindist = inf;
    for i=1:4
        dist = norm([sclon, sclat] - midp(i, :));
        if dist < mindist
            mindist = dist;
            minidx  = i;
        end
    end

    % Determine closest side
    switch minidx
        case 1, dir1 = 'north';
        case 2, dir1 = 'south';
        case 3, dir1 = 'west';
        case 4, dir1 = 'east';
    end

else
    if abs(xi - maxlon) < 1e-2
        dir1 = 'east';
    elseif abs(xi - minlon) < 1e-2
        dir1 = 'west';
    elseif abs(yi - minlat) < 1e-2
        dir1 = 'south';
    elseif abs(yi - maxlat) < 1e-2
        dir1 = 'north';
    end

    % Determine if the spacecraft is moving towards or away from the cside
    dist  = norm([sclon, sclat] - [cx, cy]);
    dist_ = norm([sclon_, sclat_] - [cx, cy]);
    if dist_ < dist
        switch dir1
            case 'north'
                dir1 = 'south';
            case 'south'
                dir1 = 'north';
            case 'east'
                dir1 = 'west';
            case 'west'
                dir1 = 'east';
        end
    end
end

% Calculate how is the spacecraft ground track position moving along the
% map. The coverage path does not only depend on the spacecraft position
% itself but also its velocity direction

switch dir1
    case {'north','south'} % Horizontal sweep
        if (sclon - sclon_) >= 0 % sc is moving leftwards
            dir2 = 'west'; % if the spacecraft is moving left (in
            % the topography map) then the coverage path should start
            % at the position furthest to the right (right -> left
            % direction)
        else % sc is moving to the right
            dir2 = 'east'; % if the spacecraft is moving right (in
            % the topography map) then the coverage path should start
            % at the position furthest to the left (left -> right
            % direction)
        end
    case {'east', 'west'}
        if (sclat - sclat_) < 0 % sc is moving upwards
            dir2 = 'north'; % if the spacecraft is moving up (in the
            % topography map) then the coverage path should start at
            % the position furthest to the bottom (down -> top
            % direction)
        else % sc is moving down
            dir2 = 'south';% if the spacecraft is moving down (in the
            % topography map) then the coverage path should start at
            % the position furthest to the top (top -> down direction)
        end
end

end