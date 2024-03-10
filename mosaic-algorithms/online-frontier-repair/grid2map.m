function map = grid2map(grid)
% Grid map: this map encloses the discretized grid. First and last rows
% and columns are NaN, as well as the empty values in the grid (points
% outside of the ROI or excluded from observation bc their associated
% allocated cell's coverage is too small)
map = cell(size(grid,1) + 2, size(grid,2) + 2);
for i=1:size(map,1)
    for j=1:size(map,2)
        if i == 1 || j == 1 || i == size(map,1) || ...
                j == size(map,2) || any(isempty(grid{i-1, j-1}))
            map{i,j} = [NaN NaN];
        else
            map{i, j} = grid{i-1, j-1};
        end
    end
end

end