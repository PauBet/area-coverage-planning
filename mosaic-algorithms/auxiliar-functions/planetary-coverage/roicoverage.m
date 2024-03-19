function [coverage, overlap] = roicoverage(target, roi, fplist)
% Provided a list of footprints on a body surface and ROIs, this function 
% computes  the cumulative coverage and overlap on each ROI, in [%]
% Overlap accounts for the % of surface covered at least more than once
% [Note:] This needs a revision, only works for a single roi!!
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         12/2023
% Revision:     1
%
% Usage:        [coverage, overlap] = roicoverage(target, roi, fplist)
%
% Inputs:
%   > target:   SPICE ID (int) or name (string) of the body
%   > roi:      roi struct that contains, at least, its boundary vertices,
%               in latitudinal coordinates (in [deg])
%   > fplist:   list of footprints (struct). See footprint for further
%               information
%
% Output:
%   > coverage: Percentage of surface that the list of provided footprints 
%               collectively cover with respect to the ROI surface,
%               in [%]
%   > overlap:  Percentage of overlap of the list of provided footprints 
%               with respect to the ROI surface, in [%]

% Pre-allocate variables
coverage = 0; % initialize coverage
overlap = 0; % initialize overlap
if isempty(fplist), return; end

% In case roi is not input as a struct but as a matrix (vertices)...
if ~isstruct(roi)
    [vertices(:, 1), vertices(:, 2)] = amsplit(roi(:, 1), roi(:, 2));
    clear roi;
    roi = struct('vertices', interppolygon(vertices));
else
    [aux(:, 1), aux(:, 2)] = amsplit(roi.vertices(:, 1), ...
        roi.vertices(:, 2));
    roi.vertices = interppolygon(aux);
end

% Cumulative coverage
for i=1:length(fplist)

    if ~isfield(roi, 'polyJ')
        roi.polyJ = polyshape();
        roi.polyO = polyshape();
    end

    % Footprint polyshape
    if isempty(fplist(i).bvertices)
        continue;
    end
    x = fplist(i).bvertices(:, 1); y = fplist(i).bvertices(:, 2);
    polyFP = polyshape(x, y);

    % Intersect footprint to cumulative coverage
    polyI  = intersect(roi.polyJ, polyFP);
    roi.polyO  = union(roi.polyO, polyI);

    % Join footprint
    roi.polyJ  = union(roi.polyJ, polyFP);
end

% Calculate ROI specific coverage
% Get total surface of the ROI
lon = roi.vertices(:, 1); lat = roi.vertices(:, 2);
RA = regionarea(target, lon, lat);
polyROI = polyshape(lon, lat); % ROI polyshape

% Intersect ROI with cumulative footprints
if ~isempty(roi.polyJ.Vertices)
    polyI = intersect(polyROI, roi.polyJ);
    if ~isempty(polyI.Vertices)
        inter = subtract(polyROI, polyI);
        coverage = (area(polyROI) - area(inter))*1e2/area(polyROI);
        %A = regionarea(target, inter.Vertices(:, 1), inter.Vertices(:, 2));
        %coverage = (RA - A)*1e2/RA;
        % lon = polyI.Vertices(:, 1); lat = polyI.Vertices(:, 2);
        % A = regionarea(target, lon, lat); % get area
        % 
        % % Compute ROI coverage
        % coverage = A*1e2/RA;
        % coverage = sum(coverage); % if there are more than one regions, the
        % % coverage of those can be simply added
        % if coverage > 100
        %     disp("pause")
        % end
    end
end

% Intersect ROI with overlap
if ~isempty(roi.polyO.Vertices)
    polyI = intersect(polyROI, roi.polyO);
    if ~isempty(polyI.Vertices)
        % This should be correct but does not work as good as the
        % above algorithm (coverage), because polyI here is very irregular,
        % hence areaint may fail to give an accurate solution
        % inter = subtract(polyROI, polyI);
        % A = regionarea(target, inter.Vertices(:, 1), inter.Vertices(:, 2));
        % overlap = (RA - A)*1e2/RA;
        
        % Alternative (not optimal either)
        inter = subtract(polyROI, polyI);
        overlap = (area(polyROI) - area(inter))*1e2/area(polyROI);
    end
end

end