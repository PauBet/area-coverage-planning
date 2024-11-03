function plotTour_m(tour, fplist, roistruct, target, varargin)

% Pre-allocate variables
filename = strcat(lower(target),'-map.jpg');
c1 = [0.72,0.27,1.00];
c2 = [0.35,0.06,0.41];
if nargin > 4
    ax = varargin{1};
    if nargin > 5 
        video = varargin{2};
    else
        video = [];
    end
else
    % Create figure
    ax = mapPlot(filename);
    video = [];
end

% Plot region(s) of interest
if ~isempty(video)
    for i=1:length(roistruct)
        roi = roistruct(i).vertices;
        poly = polyshape(roi(:, 1), roi(:, 2));
        plot(poly, 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none', ...
            'DisplayName', roistruct(i).name);
    end
end

% Check
if isempty(tour)
    return;
end

% Plot footprint list
for i=1:length(fplist)
    % Plot footprint polygon
    poly = polyshape(fplist(i).bvertices(:, 1), fplist(i).bvertices(:, 2));
    h1 = plot(poly, 'FaceColor', c1, 'FaceAlpha', 0.2, 'EdgeColor', c2, ...
        'LineWidth', .75, 'DisplayName', 'Footprint');

    % Plot coverage path
    if i > 1
        if abs(tour{i-1}(1) - tour{i}(1)) <= 180 % no coverage ...
            % path - a.m. intercept
            h2 = plot(ax, [tour{i-1}(1) tour{i}(1)], [tour{i-1}(2) tour{i}(2)],...
                'Color', 'y', 'linewidth', 1, 'DisplayName', 'Coverage path');
        end
    else
        h3 = scatter(ax, tour{i}(1), tour{i}(2), 50, 'b', "filled", '^', ...
            'DisplayName', 'Start point');
    end

    % % Plot ground track
    % t = fplist(i).t;
    % [sclon, sclat] = groundtrack(sc, t, target);
    % if i > 1
    %     h4 = scatter(ax, sclon, sclat, 8, 'c', "filled", 'DisplayName', 'Ground track');
    % else
    %     scatter(ax, sclon, sclat, 50, 'b', "filled", '^');
    % end
    
    % Save animation
    if ~isempty(video)
        drawnow;
        writeVideo(video, getframe(gcf));
    end
end

% Re-plot the ROI (for aesthetic purposes)
for i=1:length(roistruct)
    roi = roistruct(i).vertices;
    [x, y] = amsplit(roi(:, 1), roi(:, 2));
    poly = polyshape(x, y);
    h5 = plot(poly, 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceColor', 'none', ...
        'DisplayName', roistruct(i).name);
end

% Legend
legend([h1, h2, h3], 'AutoUpdate', 'off')

if ~isempty(video), close(video); end

end