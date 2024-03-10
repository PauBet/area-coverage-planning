function targetArea = topo2inst(roi, lon, lat, target, sc, inst, et)

% Pre-allocate variables
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE

% Build focal plane
[fovbounds, boresight, rotmat] = instpointing(inst, target, sc, et, lon, lat);
vertex = cspice_spkpos(sc, et, targetframe, 'NONE', target);
point = vertex + fovbounds(:, 1);
plane = cspice_nvp2pl(boresight, point);

% Intersect ROI with focal plane
spoint = zeros(size(roi, 1), 3);
for i=1:size(roi, 1)
    if ~isnan(roi(i, :))
        dir = -trgobsvec(roi(i, :), et, target, sc);
        [found, spoint(i, :)] = cspice_inrypl(vertex, dir, plane);
        if found == 0
            disp("No intersection");
        end
    else
        spoint(i, :) = nan(1, 3);
    end
end

% Transform coordinates from body-fixed to instrument frame
tArea = zeros(length(spoint), 3);
for i=1:length(spoint)
    if ~isnan(spoint(i, :))
        vpoint = -(vertex - spoint(i, :)');
        tArea(i, :) = rotmat\vpoint;
    else
        tArea(i, :) = nan(1, 3);
    end
end
targetArea = tArea(:, 1:2);

end