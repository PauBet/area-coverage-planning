function [coverage, overlap] = mapcoverage(target, fplist)
% Provided a list of footprints on a body surface, this function computes 
% the cumulative coverage and overlap, in [%]
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         12/2023
% Revision:     1
%
% Usage:        [coverage, overlap] = mapcoverage(target, fplist)
%
% Inputs:
%   > target:   SPICE ID (int) or name (string) of the body
%   > fplist:   list of footprints (struct). See footprint for further
%               information
%
% Output:
%   > coverage: Percentage of surface that the list of provided footprints 
%               collectively cover with respect to the entire body surface,
%               in [%]
%   > overlap:  Percentage of overlap of the list of provided footprints 
%               with respect to the entire body surface, in [%]

% Pre-allocate variables
polyJ  = polyshape(); % cumulative coverage
polyO  = polyshape(); % overlap
overlap = 0; % initialize overlap
coverage= 0; % initialize coverage

% Cumulative coverage
for i=1:length(fplist)

    % Footprint polyshape
    if isempty(fplist(i).bvertices)
        continue;
    end
    x = fplist(i).bvertices(:, 1); y = fplist(i).bvertices(:, 2);
    polyFP = polyshape(x, y);

    % Intersect footprint to cumulative coverage
    polyI  = intersect(polyJ, polyFP);
    polyO  = union(polyO, polyI);

    % Join footprint
    polyJ  = union(polyJ, polyFP);
end

% Coverage of all footprints
if ~isempty(polyJ.Vertices)
    lon = polyJ.Vertices(:, 1);
    lat = polyJ.Vertices(:, 2);
    [~, coverage] = regionarea(target, lon, lat);
end
coverage = sum(coverage); % if there are more than one regions, the 
% coverage of those can be simply added

% Coverage of overlap
if ~isempty(polyO.Vertices)
    lon = polyO.Vertices(:, 1);
    lat = polyO.Vertices(:, 2);
    [~, overlap] = regionarea(target, lon, lat);
end
overlap = sum(overlap);

end