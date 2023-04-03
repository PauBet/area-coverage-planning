% Define tile sizes
tile_sizes = [1, 2, 1, 3, 1, 1, 2, 1];

% Calculate the cumulative sum of tile sizes
cumulative_sizes = cumsum(tile_sizes);

% Calculate the total size of the grid
grid_size = cumulative_sizes(end);

% Calculate the optimal stare points
stare_points = zeros(1, length(tile_sizes));
for i = 1:length(tile_sizes)
    if i == 1
        stare_points(i) = tile_sizes(i) / 2;
    else
        stare_points(i) = (cumulative_sizes(i-1) + tile_sizes(i) / 2) ...
            - ((grid_size - cumulative_sizes(end)) / (length(tile_sizes) - 1)) * (i - 1);
    end
end

% Plot the tile grid and stare points
figure
hold on
for i = 1:length(tile_sizes)
    rectangle('Position', [i-tile_sizes(i)/2, 0, tile_sizes(i), 1], 'FaceColor', 'k')
end
plot(stare_points, zeros(1, length(tile_sizes)), 'ro')
xlim([0.5, length(tile_sizes)+0.5])
axis equal