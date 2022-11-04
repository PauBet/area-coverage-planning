clear; clc; close all;

addpath(genpath(pwd));

%%
startTime = cspice_str2et('2011 SEP 30 2:00:00.000 TDB');
endTime   = cspice_str2et('2011 SEP 30 2:30:00.000 TDB');
step = 60;
instName = 'DAWN_FC2';
scName = 'DAWN';
targetName = 'VESTA';
method = 'ELLIPSOID';
[~, targetFrame, ~] = cspice_cnmfrm(targetName); % body-fixed frame

%%

% polygon (from the frontier vertices)
areaPoints = [ -120 30;
               -115 25;
               -110 20;
               -110 15;
               -110 10;
               -115 10;
               -120 10];

areaPoints = [10 -10;
              10 -50;
              60 -60;
              60 -10;];

poly1 = polyshape(areaPoints(:,1), areaPoints(:,2));

%%
ax2 = mapPlot('vesta-map.png');
fig2 = gcf;
set(gcf,'units','normalized','OuterPosition',[0.4307,0.3144,0.5656,0.6565]);
plot(ax2, polyshape(areaPoints), 'FaceColor', [0.93,0.69,0.13])

% discretize the target area
[gamma(1), gamma(2)] = centroid(poly1);
fprint0 = footprint(gamma(1), gamma(2), instName, targetName, scName, ...
    startTime, method, [], [], []); % compute the observation's footprint
olapx = 0;
olapy = 0;
[grid, vlon, vlat] = grid2D(fprint0.sizex, fprint0.sizey, olapx, ...
    olapy, gamma, areaPoints);

% build the map
map = cell(size(grid,1) + 2, size(grid,2) + 2);
for i=1:size(map,1)
    for j=1:size(map,2)
        if i == 1 || j == 1 || i == size(map,1) || j == size(map,2)
            map{i,j} = [NaN, NaN];
        else
            map{i, j} = grid{i-1, j-1};
        end
    end
end

% Initial frontier vertices
F = {};
for i=1:size(map,1)
    for j=1:size(map, 2)
        n = getNeighbors(i, j, map);
        if length(n) < 8 && ~isempty(n)
            F{end + 1} = [map{i,j}];
        end
    end
end

% tour - rectangular Boustrophedon decomposition
tour = {};
bearing = true;
for i=1:size(map,1)
    for j=1:size(map,2)
        if bearing
            if ~isnan(map{i,j})
                tour{end + 1} = map{i,j};
            end
        else
            if ~isnan(map{i,size(map,2) + 1 - j})
                tour{end + 1} = map{i,size(map,2) + 1 - j};
            end
        end
        bearing = not(bearing);
    end
end

A = {};
t = startTime;
while ~isempty(poly1.Vertices) && t <= endTime
    [N, X, F] = updateGrid(areaPoints, tour, F, map); % seed fill algorithm
    for i=1:length(X)
        for j=1:length(tour)
            if isequal(X{i}, tour{j})
                tour(j) = [];
            end
        end
    end
    for i=1:length(N)
        if ~isequal(N{i}, a)
            tour{end+1} = N{i};
        end
    end
    % Compute the footprint of each point in the tour successively and
    % subtract the corresponding area from the target polygon
    a = tour{1}; % observation
    tour(1) = []; % delete this observation from the planned tour
    A{end + 1} = a; % add it in the list of planned observations
    fprinti = footprint(a(1), a(2), instName, targetName, scName, ...
        t, method, [], [], []); % compute the observation's footprint
    if ~isempty(fprinti.bbox)
        poly2 = polyshape(fprinti.bbox); % create footprint polygon
        poly1 = subtract(poly1, poly2); % update uncovered area
        vertices = poly1.Vertices;
    end
    t = t + step;

    %% Plots

    % Footprint plot in Figure 2
    plot(ax2, poly2, 'FaceColor', [0.00,0.45,0.74])
    if length(A) > 1
        plot(ax2, [A{end-1}(1) A{end}(1)], ...
            [A{end-1}(2) A{end}(2)], 'w-', 'linewidth', 1)
    end
    drawnow

    % Spacecraft ground track position in Figure 2
    sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', targetName,...
        t, targetFrame, 'NONE', scName);
    [~, sclon, sclat] = cspice_reclat(sctrack);
    plot(sclon*cspice_dpr, sclat*cspice_dpr, 'y^',...
        'MarkerSize', 6)
end