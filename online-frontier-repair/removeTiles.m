function [tour, map] = removeTiles(tour, map, tiles)
% This function removes disposable observation points in a planned tour
% observations (and the grid discretization)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        [tour, map] = removeTiles(tour, map, tiles)
%
% Inputs:
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > tiles:        cell matrix of the disposable observation points to be
%                   removed from 'tour' and 'map'
% 
% Outputs:
%   > tour:        updated cell matrix of the successive planned
%                  observations
%   > map:         updated cell matrix of grid points

indel = [];
for i=1:length(tiles)
    % For each observation point in the removal list...

    for j=1:length(tour)
        if norm(tiles{i} - tour{j}) < 1e-5
            % Find index position in 'tour'
            indel = [indel j];
        end
    end
    for ii=1:size(map, 1)
        for jj=1:size(map, 2)
            if norm(map{ii, jj} - tiles{i}) < 1e-5

                % Remove elements from the grid ([NaN NaN])
                map{ii, jj} = [NaN NaN];
            end
        end
    end
end

% Remove elements from tour
tour(indel) = [];

end