function varargout = grid2D(w, h, ovlapx, ovlapy, gamma, targetArea)
 %% Paula
gridPoints = [];
figure
plot(polyshape(targetArea(:,1), targetArea(:,2)))
hold on;
plot(gamma(1), gamma(2), 'b*')
gridPoints = floodFillAlgorithm(w, h, ovlapx, ovlapy, gamma, targetArea, gridPoints, '4fill');
xlim([min(targetArea(:,1)) max(targetArea(:,1))])
ylim([min(targetArea(:,2)) max(targetArea(:,2))])

% The elements of gridPoints are sorted by latitude (+ to -)
sortedGrid = sortrows(gridPoints, -2);
% Build a matrix that will sort the gridPoints elements by latitude and
% longitude according to the following structure:
%
%               longitude
%               (-) --------> (+)                   
%  latitude (+) [a11]  [a12] ⋯
%            ¦  [a21]
%            ¦    ⋮
%            ∨
%           (-)             
%             
uniqueLat = unique(sortedGrid(:,2)); % get the different latitude values
% unique check
indlat = find(abs(diff(uniqueLat)) < 1e-5);
uniqueLat(indlat) = [];
uniqueLon = unique(sortedGrid(:,1)); % get the different longitude values
% unique check
indlat = find(abs(diff(uniqueLon)) < 1e-5);
uniqueLon(indlat) = [];

matrixGrid = cell(length(uniqueLat), length(uniqueLon));
for i=1:size(matrixGrid,1)
    for j=1:size(matrixGrid,2)
        matrixGrid{i,j} = [NaN, NaN];
    end
end

%
for i=1:length(uniqueLat)
    lat = uniqueLat(length(uniqueLat) + 1 - i);
    indlat = (abs(sortedGrid(:,2) - lat) < 1e-5);
    mrow = sort(sortedGrid(indlat));
    for j=1:length(mrow)
        indlon = find(abs(uniqueLon - mrow(j)) < 1e-5);
        matrixGrid{i, indlon} = [mrow(j), lat];
    end
end

%%
varargout{1} = matrixGrid;
if nargout > 1
    varargout{2} = uniqueLon;
    varargout{3} = uniqueLat;
end

end