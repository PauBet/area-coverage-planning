function cPoint = closestTargetCorner(targetName, scName, t, map)
% [Description]
%%
[~, targetFrame, ~] = cspice_cnmfrm(targetName);
% Compute spacecraft sub-observer point to see what side of the target
% area is closer to it
scSubObs = cspice_subpnt('NEAR POINT/ELLIPSOID', targetName, t,...
    targetFrame, 'NONE', scName);
[~, sclon, sclat] = cspice_reclat(scSubObs);
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr;
% Compute the closest target-area corner to the spacecraft's
% sub-observer point
rows = [1 size(map,1)];
cols = [1 size(map,2)];
mindist = inf;
for i=1:length(rows)
    for j=1:length(rows)
        dist = norm([sclon, sclat] - map{rows(i),cols(j)});
        if dist < mindist
            indrow = rows(i); indcol = cols(j);
            mindist = dist;
        end
    end
end
cPoint = map{indrow, indcol};
end