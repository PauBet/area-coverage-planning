function outputData = topo2inst(inputdata, lon, lat, target, sc, inst, et)
% This function transforms a set of points from the topographic coordinate 
% system(latitude and longitude on the target body) to the instrument frame 
% coordinates.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        outputData = topo2inst(inputdata, lon, lat, target, sc, inst, et)
%
% Inputs:
%   > inputdata:    cell array or matrix of points in topographic 
%                   coordinates to be transformed. Each point is a row with 
%                   [longitude, latitude] format
%   > lon:          longitude of the observation point or area center, in
%                   [deg]
%   > lat:          latitude of the observation point or area center, in
%                   [deg]
%   > target:       string name of the target body
%   > sc:           string name of the spacecraft
%   > inst:         string name of the instrument
%   > et:           ephemeris time, TDB seconds past J2000 epoch
% 
% Outputs:
%   > outputData:   A cell array or matrix of the input points transformed 
%                   to the instrument frame coordinates. The format of the 
%                   output matches the input (cell array or matrix)

% Handle input data in cell format, ensuring all empty entries are replaced
% with [NaN, NaN]
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
% Create a plane based on the boresight and a point in the focal plane
plane = cspice_nvp2pl(boresight, point);

% For each topographic point, find its intersection with the focal plane
spoint = zeros(size(topoPoints, 1), 3);
for i=1:size(topoPoints, 1)
    if ~isnan(topoPoints(i, :))
        dir = -trgobsvec(topoPoints(i, :), et, target, sc);
        [found, spoint(i, :)] = cspice_inrypl(vertex, dir, plane);
        if found
            emnang = emissionang(topoPoints(i, :), et, target, sc);
            if emnang >= 90, found = 0; end
        end
        if found == 0
            spoint(i, :) = nan(1, 3);
            %disp("No intersection");
        end
    else
        spoint(i, :) = nan(1, 3);
    end
end

% Transform coordinates from body-fixed to instrument frame
tArea = zeros(size(spoint, 1), 3);
for i=1:size(spoint, 1)
    if ~isnan(spoint(i, :))
        vpoint = -(vertex - spoint(i, :)'); % vector from spacecraft to 
        % intersection point
        tArea(i, :) = rotmat\vpoint; % apply inverse rotiation to transform
        % to instrument frame
    else
        tArea(i, :) = nan(1, 3);
    end
end
instcoord = tArea(:, 1:2); % extract 2D instrument frame coordinates

% Prepare output data matching the format of the input, i.e., cell array or
% matrix
outputData = cell(size(inputdata));
if iscell(inputdata)
    for k=1:length(ii)
        if ~isnan(instcoord(k))
            outputData{ii(k), jj(k)} = instcoord(k, :);
        end
    end
else
    instcoord(isnan(instcoord(: ,1)), :) = [];
    [aux(:, 1), aux(:, 2)] = sortcw(instcoord(:, 1), instcoord(:, 2));
    outputData = aux;
end

end