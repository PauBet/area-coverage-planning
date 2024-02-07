function outputData = topo2inst(inputdata, lon, lat, target, sc, inst, et)

if iscell(inputdata)
    for i=1:size(inputdata, 1)
        for j=1:size(inputdata, 2)
            if ~isempty(inputdata{i, j})
                aux{i, j} = inputdata{i, j};
            else
                aux{i, j} = [NaN NaN];
            end
        end
    end
    topoPoints = vertcat(aux{:});
    [ii, jj] = ind2sub(size(inputdata), 1:numel(inputdata));
else
    topoPoints = inputdata;
end

% Pre-allocate variables
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE

% Build focal plane
[fovbounds, boresight, rotmat] = instpointing(inst, target, sc, et, lon, lat);
vertex = cspice_spkpos(sc, et, targetframe, 'NONE', target);
point = vertex + fovbounds(:, 1);
plane = cspice_nvp2pl(boresight, point);

% Intersect topoPoints with focal plane
spoint = zeros(size(topoPoints, 1), 3);
for i=1:size(topoPoints, 1)
    if ~isnan(topoPoints(i, :))
        dir = -trgobsvec(topoPoints(i, :), et, target, sc);
        [found, spoint(i, :)] = cspice_inrypl(vertex, dir, plane);
        if found == 0
            disp("No intersection");
        end
    else
        spoint(i, :) = nan(1, 3);
    end
end

% Transform coordinates from body-fixed to instrument frame
tArea = zeros(size(spoint, 1), 3);
for i=1:size(spoint, 1)
    if ~isnan(spoint(i, :))
        vpoint = -(vertex - spoint(i, :)');
        tArea(i, :) = rotmat\vpoint;
    else
        tArea(i, :) = nan(1, 3);
    end
end
instcoord = tArea(:, 1:2);

%
outputData = cell(size(inputdata));
if iscell(inputdata)
    for k=1:length(ii)
        if ~isnan(instcoord(k))
            outputData{ii(k), jj(k)} = instcoord(k, :);
        end
    end
else
    outputData = instcoord;
end

end