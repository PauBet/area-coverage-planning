clear; clc; close all;

addpath(genpath(pwd));

%%
inputkernels;
addpath('/Users/paulabetriu/Desktop/GitHub/RESSlib');
initSPICEv(fullK(METAKR));

%%
startTime = cspice_str2et('2011 SEP 30 1:48:00.000 TDB');
endTime   = cspice_str2et('2011 SEP 30 2:00:00.000 TDB');
step = 60;
inst = 'DAWN_FC2';
sc = 'DAWN';
target = 'VESTA';
method = 'ELLIPSOID';
[~, targetFrame, ~] = cspice_cnmfrm(target); % body-fixed frame

%
t = startTime;

%%

% polygon (from the frontier vertices)
roi = [ -120 30;
        -115 25;
        -110 20;
        -110 15;
        -110 10;
        -115 10;
        -120 10];

roi = [10 -10;
       10 -50;
       60 -60;
       60 -10;];

poly1 = polyshape(roi(:,1), roi(:,2));

%%
ax2 = mapPlot('vesta-map.png');
fig2 = gcf;
set(gcf,'units','normalized','OuterPosition',[0.4307,0.3144,0.5656,0.6565]);
plot(ax2, polyshape(roi), 'FaceColor', [0.93,0.69,0.13])

% discretize the target area
[gamma(1), gamma(2)] = centroid(poly1);
fprint0 = footprint(gamma(1), gamma(2), startTime, inst, ...
    sc, target, 0); % compute the observation's footprint
olapx = 0;
olapy = 0;
grid = grid2D(fprint0, olapx, olapy, gamma, roi);

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

% Closest polygon side to the spacecraft's ground track position (this
% will determine the coverage path in planSidewinderTour)
cside = closestSide(target, sc, t, roi);

% Sorted list of grid points according to the sweeping/coverage path
% (see Boustrophedon decomposition)
tour = planSidewinderTour(target, sc, t, roi, fprint0, gamma,...
    olapx, olapy, cside);

A = {};
t = startTime;
exit = false;
while t <= endTime
    [N, X, F] = updateGrid(roi, tour, F, map); % seed fill algorithm
    indel = [];
    for i=1:length(X)
        for j=1:length(tour)
            if isequal(X{i}, tour{j})
                indel = [indel j];
            end
        end
    end
    tour(indel) = [];
    if ~isempty(tour)
    % Compute the footprint of each point in the tour successively and
    % subtract the corresponding area from the target polygon
    a = tour{1}; % observation
    for i=1:length(N)
        if ~isequal(N{i}, a)
            tour{end+1} = N{i};
        end
    end
    tour(1) = []; % delete this observation from the planned tour
    A{end + 1} = a; % add it in the list of planned observations
    fprinti = footprint(a(1), a(2), t, inst, ...
    sc, target, 0); % compute the observation's footprint
    if ~isempty(fprinti.bvertices)
        poly2 = polyshape(fprinti.bvertices); % create footprint polygon
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
    sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
        t, targetFrame, 'NONE', sc);
    [~, sclon, sclat] = cspice_reclat(sctrack);
    plot(sclon*cspice_dpr, sclat*cspice_dpr, 'y^',...
        'MarkerSize', 6)

    roi = poly1.Vertices;

    else
        break;
    end
end