function matrixGrid = sortLatLonGrid(gridPoints)

% Pre-allocate variables
matrixGrid = {};

if ~isempty(gridPoints)

    % Sort grid points
    sortedGrid = sortrows(gridPoints, -2); % the elements of gridPoints are
    % sorted by latitude (+ to -)
    uniqueLat = unique(sortedGrid(:,2)); % get the different latitude values
    ind = abs(diff(uniqueLat)) < 1e-5; % double check that there are
    % no "similar" latitude values (it may happen)
    uniqueLat(ind) = [];
    uniqueLon = unique(sortedGrid(:,1)); % get the different longitude values
    % unique check
    ind = abs(diff(uniqueLon)) < 1e-5; % double check that there are
    % no "similar" longitude values (it may happen)
    uniqueLon(ind) = [];

    % Sort and rotate the grid points and insert them in the grid matrix
    matrixGrid = cell(length(uniqueLat), length(uniqueLon));
    for i=1:length(uniqueLat)
        % We will sweep across the grid by, first, latitude and, second,
        % longitude
        lat = uniqueLat(length(uniqueLat) + 1 - i);
        indlat = (abs(sortedGrid(:,2) - lat) < 1e-5);
        mrow = sort(sortedGrid(indlat));
        for j=1:length(mrow)
            indlon = abs(uniqueLon - mrow(j)) < 1e-5;
            lon = mrow(j);
            matrixGrid{i, indlon} = [cx, cy]' + transpose(rotmat)*...
                ([lon lat]' - [cx, cy]'); % rotate values to the original roi
            % orientation
        end
    end
end
end