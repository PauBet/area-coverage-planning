function map = removeTiles(map, tiles)
% This function removes disposable observation points within the grid
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         06/2023
% 
% Usage:        map = removeTiles(map, tiles)
%
% Inputs:
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > tiles:        cell matrix of the disposable observation points to be
%                   removed from 'tour' and 'map'
% 
% Outputs:
%   > map:         updated cell matrix of grid points

for i=1:length(tiles)
    % For each observation point in the removal list...

    for ii=1:size(map, 1)
        for jj=1:size(map, 2)
            if norm(map{ii, jj} - tiles{i}) < 1e-5

                % Remove elements from the grid ([NaN NaN])
                map{ii, jj} = [NaN NaN];
            end
        end
    end
end

end