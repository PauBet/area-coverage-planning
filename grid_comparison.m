clc; close all; clear all;
% Revision of Sidewinder:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see
% Re-implementation of sidewinder function with these new feature (and
% other code improvements)

% Load mission info (kernels, SPICE ids, etc.)
input_data;
%inittime = cspice_str2et('1998 MAR 29 12:53:00.000 TDB'); % closest approach
inittime = cspice_str2et('1998 MAR 29 12:10:00.000 TDB'); % closest approach
stoptime = inittime + 24*3600;
tcadence = 15; % observation cadence [s]
olapx = 20; olapy = 20; % overlap in [%]
speedUp = 0;

% Define ROI (Annwn Regio)
roi = [25 45;
       25 15;
       -5  15;
       -5  45;];

roi = [25  25;
       25  -5;
       -5  -5;
       -5  25;];

% Point camera at ROI's centroid
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(1).vertices = roi;
roistruct(1).cpoint = [cx, cy];
roistruct(1).name = "Region of Interest";

%% Topographical plane grid
% Initial 2D grid layout discretization: the instrument's FOV is going
% to be projected onto the uncovered area's centroid and the resulting
% footprint shape is used to set the grid spatial resolution
[gamma(1), gamma(2)] = centroid(polyshape(roi(:,1),roi(:,2)));
fpref = footprint(inittime, inst, sc, target, 'highres', ...
            gamma(1), gamma(2)); % centroid footprint
grid = grid2D(fpref, olapx, olapy, gamma, roi);

% Plot
% Pre-allocate variables
figure
hold on; box on; grid minor; axis equal;
c1 = [0.93,0.69,0.13];
c2 = [0,0,0];
% Create figure
%ax = mapPlot(filename);
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        if ~isempty(grid{i, j})
            lon = grid{i, j}(1); lat = grid{i, j}(2);
            fprinti = fpref;
            fprinti.bvertices(:, 1) = fpref.bvertices(:, 1) + lon - cx;
            fprinti.bvertices(:, 2) = fpref.bvertices(:, 2) + lat - cy;
            h1 = plot(polyshape(fprinti.bvertices(:, 1), ...
                fprinti.bvertices(:, 2)), 'FaceColor', c1, 'FaceAlpha', 0.2, ...
                'EdgeColor', c2, 'LineWidth', .75, 'DisplayName', 'Ref. Tile');
            h2 = plot(lon, lat, '^', 'MarkerEdgeColor', [0.66,0.2,0.23], ...
                'MarkerSize', 5, 'DisplayName', 'Stare point', 'LineWidth', .75);
            hold on;
        end
    end
end
for i=1:length(roistruct)
    roi = roistruct(i).vertices;
    poly = polyshape(roi(:, 1), roi(:, 2));
    h3 = plot(poly, 'EdgeColor', [0.66,0.2,0.23], 'LineWidth', 1.5, 'FaceColor', 'none', ...
        'DisplayName', roistruct(i).name);
end

% Legend
legend([h1, h2, h3], 'AutoUpdate', 'off', 'Location', 'northwest', 'NumColumns', 3)

% % Zoom in
xlim([-30 50])
ylim([-15 35])
xtick = -30:20:50;
ytick = -15:10:35;

% x tick label
xtickstr = string();
for i=1:length(xtick)
    if xtick(i) < 0 && xtick(i) > -180
        xtickstr(i) = strcat(num2str(-xtick(i)), 'º', 'W');
    elseif xtick(i) > 0 && xtick(i) < 180
        xtickstr(i) = strcat(num2str(xtick(i)), 'º', 'E');
    else
        xtickstr(i) = strcat(num2str(abs(xtick(i))), 'º');
    end
end

% y tick label
ytickstr = string();
for i=1:length(ytick)
    if ytick(i) < 0
        ytickstr(i) = strcat(num2str(-ytick(i)), 'º', 'S');
    elseif ytick(i) > 0
        ytickstr(i) = strcat(num2str(ytick(i)), 'º', 'N');
    else
        ytickstr(i) = strcat(num2str(ytick(i)), 'º');
    end
end
set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
    'YTickLabel', ytickstr, ...
    'Layer', 'top', 'GridColor', 'k', 'GridAlpha', 1, 'GridLineStyle', ':', ...
    'FontSize', 20, 'InnerPosition', [0.1378, 0.15, 0.7672, 0.807]);
xlabel('Planetocentric longitude [º]')
ylabel('Planetocentric latitude [º]')

% Save figure [PDF]
figpath = '.';
set(gcf, 'Units', 'inches', 'Position', [3,3,8,5.4]);
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
filename = fullfile(figpath, 'floodfill');
print(gcf, filename, '-dpdf', '-r600')

% Plot
% Pre-allocate variables
figure
hold on; box on; grid minor; axis equal;
c1 = [0.93,0.69,0.13];
c2 = [0,0,0];
% Create figure
%ax = mapPlot(filename);
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        if ~isempty(grid{i, j})
            lon = grid{i, j}(1); lat = grid{i, j}(2);
            fprinti = footprint(inittime, inst, sc, target, 'highres', lon, lat);
            h1 = plot(polyshape(fprinti.bvertices(:, 1), ...
                fprinti.bvertices(:, 2)), 'FaceColor', c1, 'FaceAlpha', 0.2, ...
                'EdgeColor', c2, 'LineWidth', .75, 'DisplayName', 'Footprint');
            h2 = plot(lon, lat, '^', 'MarkerEdgeColor', [0.66,0.2,0.23], ...
                'MarkerSize', 5, 'DisplayName', 'Stare point', 'LineWidth', .75);
            hold on;
        end
    end
end
for i=1:length(roistruct)
    roi = roistruct(i).vertices;
    poly = polyshape(roi(:, 1), roi(:, 2));
    h3 = plot(poly, 'EdgeColor', [0.66,0.2,0.23], 'LineWidth', 1.5, 'FaceColor', 'none', ...
        'DisplayName', roistruct(i).name);
end

% Legend
legend([h1, h2, h3], 'AutoUpdate', 'off', 'Location', 'northwest', 'NumColumns', 3)

% Zoom in
xlim([-30 50])
ylim([-15 35])
xtick = -30:20:50;
ytick = -15:10:35;

% x tick label
xtickstr = string();
for i=1:length(xtick)
    if xtick(i) < 0 && xtick(i) > -180
        xtickstr(i) = strcat(num2str(-xtick(i)), 'º', 'W');
    elseif xtick(i) > 0 && xtick(i) < 180
        xtickstr(i) = strcat(num2str(xtick(i)), 'º', 'E');
    else
        xtickstr(i) = strcat(num2str(abs(xtick(i))), 'º');
    end
end

% y tick label
ytickstr = string();
for i=1:length(ytick)
    if ytick(i) < 0
        ytickstr(i) = strcat(num2str(-ytick(i)), 'º', 'S');
    elseif ytick(i) > 0
        ytickstr(i) = strcat(num2str(ytick(i)), 'º', 'N');
    else
        ytickstr(i) = strcat(num2str(ytick(i)), 'º');
    end
end
set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
    'YTickLabel', ytickstr, ...
    'Layer', 'top', 'GridColor', 'k', 'GridAlpha', 1, 'GridLineStyle', ':', ...
    'FontSize', 20, 'InnerPosition', [0.1378, 0.15, 0.7672, 0.807]);
xlabel('Planetocentric longitude [º]')
ylabel('Planetocentric latitude [º]')

% Save figure [PDF]
figpath = '.';
set(gcf, 'Units', 'inches', 'Position', [3,3,8,5.4]);
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
filename = fullfile(figpath, 'grid_old');
print(gcf, filename, '-dpdf', '-r600')

%% Focal plane grid
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
auxroi = interppolygon(roi);
targetArea = topo2inst(auxroi, cx, cy, target, sc, inst, inittime);
% Build reference tile
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
[angle, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
fpref.width = width;
fpref.height = height;
fpref.angle = angle;
grid = grid2D(fpref, olapx, olapy, [0, 0], targetArea);

% Transform and project coordinates
grid = inst2topo(grid, cx, cy, target, sc, inst, inittime);

% Plot
% Pre-allocate variables
filename = strcat(lower(target),'-map.jpg');
figure
hold on; box on; grid minor; axis equal;
c1 = [0.93,0.69,0.13];
c2 = [0,0,0];
% Create figure
%ax = mapPlot(filename);
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        if ~isempty(grid{i, j})
            lon = grid{i, j}(1); lat = grid{i, j}(2);
            fprinti = footprint(inittime, inst, sc, target, 'highres', lon, lat);
            h1 = plot(polyshape(fprinti.bvertices(:, 1), ...
                fprinti.bvertices(:, 2)), 'FaceColor', c1, 'FaceAlpha', 0.2, ...
                'EdgeColor', c2, 'LineWidth', .75, 'DisplayName', 'Footprint');
            h2 = plot(lon, lat, '^', 'MarkerEdgeColor', [0.66,0.2,0.23], ...
                'MarkerSize', 5, 'DisplayName', 'Stare point', 'LineWidth', .75);
            hold on;
        end
    end
end
for i=1:length(roistruct)
    roi = roistruct(i).vertices;
    poly = polyshape(roi(:, 1), roi(:, 2));
    h3 = plot(poly, 'EdgeColor', [0.66,0.2,0.23], 'LineWidth', 1.5, 'FaceColor', 'none', ...
        'DisplayName', roistruct(i).name);
end

% Legend
legend([h1, h2, h3], 'AutoUpdate', 'off', 'Location', 'northwest', 'NumColumns', 3)

% % Zoom in
xlim([-30 50])
ylim([-15 35])
xtick = -30:20:50;
ytick = -15:10:35;

% x tick label
xtickstr = string();
for i=1:length(xtick)
    if xtick(i) < 0 && xtick(i) > -180
        xtickstr(i) = strcat(num2str(-xtick(i)), 'º', 'W');
    elseif xtick(i) > 0 && xtick(i) < 180
        xtickstr(i) = strcat(num2str(xtick(i)), 'º', 'E');
    else
        xtickstr(i) = strcat(num2str(abs(xtick(i))), 'º');
    end
end

% y tick label
ytickstr = string();
for i=1:length(ytick)
    if ytick(i) < 0
        ytickstr(i) = strcat(num2str(-ytick(i)), 'º', 'S');
    elseif ytick(i) > 0
        ytickstr(i) = strcat(num2str(ytick(i)), 'º', 'N');
    else
        ytickstr(i) = strcat(num2str(ytick(i)), 'º');
    end
end
set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
    'YTickLabel', ytickstr, ...
    'Layer', 'top', 'GridColor', 'k', 'GridAlpha', 1, 'GridLineStyle', ':', ...
    'FontSize', 20);
xlabel('Planetocentric longitude [º]')
ylabel('Planetocentric latitude [º]')

% Save figure [PDF]
figpath = '.';
set(gcf, 'Units', 'inches', 'Position', [3,3,8,5.4]);
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
filename = fullfile(figpath, 'grid_new');
print(gcf, filename, '-dpdf', '-r600')