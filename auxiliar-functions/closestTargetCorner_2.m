function cPoint = closestTargetCorner(target, sc, t, map)
% [Description]
%%
[~, targetFrame, ~] = cspice_cnmfrm(target);
% Compute spacecraft sub-observer point to see what side of the target
% area is closer to it
scSubObs = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t,...
    targetFrame, 'NONE', sc);
[~, sclon, sclat] = cspice_reclat(scSubObs);
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr;
% Compute the closest target-area corner to the spacecraft's
% sub-observer point
mindist = inf;
for i=1:size(map, 1)
        if ~isempty(map(i, 1))
            dist = norm([sclon, sclat] - [map(i,1), map(i,2)]);
            if dist < mindist
                indrow = i; 
                mindist = dist;
            end
        end
end
cPoint = map(indrow, :);
end