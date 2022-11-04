function closestSide = getClosestSide(targetName, scName, t, vertices)

[~, targetFrame, ~] = cspice_cnmfrm(targetName);
% Compute spacecraft sub-observer point to see what side of the target
% area is closer to it
scSubObs = cspice_subpnt('NEAR POINT/ELLIPSOID', targetName, t,...
    targetFrame, 'NONE', scName);
[~, sclon, sclat] = cspice_reclat(scSubObs);
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr;
% Boundary box
bbox = computeBoundingBox(vertices);
% Compute the closest side of the boundary box to the spacecraft
% sub-observer point
minlondist = 360;
minlatdist = 180;
for i=1:length(vertices)
    londist = sclon - vertices(i,1);
    latdist = abs(sclat - vertices(i,2));

    if abs(londist) < minlondist
        minlondist = londist;
        if londist > 0
            if (sclon - bbox.maxlon) > 0
                closestSide = 'west';
            else
                closestSide = 'east';
            end
        elseif londist < 0
            if (sclon - bbox.minlon) > 0
                closestSide = 'east';
            else
                closestSide = 'west';
            end
        end
    end

    if abs(latdist) < minlatdist
        minlatdist = latdist;
        if latdist > 0
            if (sclat - bbox.maxlat) > 0
                closestSide = 'north';
            else
                closestSide = 'south';
            end
        elseif latdist < 0
            if (sclat - bbox.minlat) > 0
                closestSide = 'south';
            else
                closestSide = 'north';
            end
        end
    end
end
end