function grid = map2grid(map)
grid = map(2:end-1, 2:end-1);
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        if any(isnan(grid{i, j}), 'all')
            grid{i, j} = [];
        end
    end
end
end