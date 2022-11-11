function cside = closestSide(target, sc, t, roi)
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
%
% Parameters
[~, targetFrame, ~] = cspice_cnmfrm(target); % target body-fixed reference 
% frame

% Compute spacecraft sub-observer point to see what side of the target
% area is closer to it
% Assumption: Tri-axial ellipsoid to model the target surface
subobs = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t,...
    targetFrame, 'NONE', sc);
[~, sclon, sclat] = cspice_reclat(subobs); % latitudinal coordinates
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr; % [rad] to [deg]

% Compute the closest side of the boundary box to the spacecraft 
% sub-observer point
bbox = smallestBoundingBox(roi(:, 1), roi(:, 2)); % get roi's bounding box
minlondist = 360;
minlatdist = 180;
for i=1:length(roi)
    londist = sclon - roi(i,1);
    latdist = abs(sclat - roi(i,2));

    if abs(londist) < minlondist
        minlondist = londist;
        if londist > 0
            if (sclon - bbox.maxlon) > 0
                cside = 'left';
            else
                cside = 'right';
            end
        elseif londist < 0
            if (sclon - bbox.minlon) > 0
                cside = 'right';
            else
                cside = 'left';
            end
        end
    end

    if abs(latdist) < minlatdist
        minlatdist = latdist;
        if latdist > 0
            if (sclat - bbox.maxlat) > 0
                cside = 'up';
            else
                cside = 'down';
            end
        elseif latdist < 0
            if (sclat - bbox.minlat) > 0
                cside = 'down';
            else
                cside = 'up';
            end
        end
    end
end
end