function A = gridNibbler(startTime, endTime, tobs, inst, sc, target, ...
    roi, olapx, olapy, videosave)
%% [Description]

%% Initialization of variables
A = {};
scPos = [];
x = roi(:,1); y = roi(:,2);
poly1 = polyshape(x,y);
[~, targetFrame, ~] = cspice_cnmfrm(target);
theta = 0; % temporal!!!!!

%% Figure 2.
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex objects
ax2 = mapPlot('vesta-map.png');
fig2 = gcf;
set(gcf,'units','normalized','OuterPosition',[0.4307,0.3144,0.5656,0.6565]);
plot(ax2, polyshape(roi), 'FaceColor', [0.93,0.69,0.13])

% Animation of figure 2
v2 = VideoWriter('topography_map_sidewinder');
v2.FrameRate = 2;
open(v2);
writeVideo(v2,getframe(fig2));

%% First iteration
t = startTime;
scPos(:,end+1) = cspice_spkpos(sc, t, targetFrame, 'NONE', target);
% Assumption: we choose to start in one of the target's area corners,
% the one closest to the spacecraft's ground track position
[gamma(1), gamma(2)] = centroid(polyshape(roi(:,1),roi(:,2)));
fp0 = footprint(gamma(1), gamma(2), t, inst,  sc, target, 0); 
map = grid2D(fp0, olapx, olapy, gamma, roi);
best = closestTargetCorner(target, sc, t, map);

%% Local Grid Nibbler
fprinti = fp0; % initial footprint
poly2 = polyshape(fprinti.bvertices); % footprint polygon
while ~iszero(roi) && t <= endTime
    
    inter = subtract(poly1, poly2);
    areaI = area(inter);
    areaT = area(poly1);
    areaInter = areaT - areaI;
    if areaInter/area(poly2) < 0.2
        newStart = closestTargetCorner(bestOld, t);
        bestNew = checkNeighbors(newStart, t);
        if score(bestNew, center) > score(best, center)
            best = newStart;
        end
    end
    
    % Footprint projection and subtraction from the ROI
    fprinti = footprint(best(1), best(2), t, inst, sc, target,...
        theta); % compute the observation's footprint
    poly2 = polyshape(fprinti.bvertices); % create footprint polygon
    poly1 = subtract(poly1, poly2); % update uncovered area
    roi = poly1.Vertices;

    %% Plots
    % Footprint plot in Figure 2
    plot(ax2, poly2, 'FaceColor', [0.00,0.45,0.74])
    if length(A) > 1
        plot(ax2, [A{end-1}(1) A{end}(1)], ...
            [A{end-1}(2) A{end}(2)], 'w-', 'linewidth', 1)
    end
    drawnow

    % Figure showing the footprint projection onto the body surface
    footprint3Dprojection(fprinti, videosave)

    % Starting from the previous tour point, all of its 8 neighbours are
    % going to be checked (and scored)
    bestOld = best;
    [center(1), center(2)] = centroid(polyshape(roi(:,1), roi(:,2)));
    best = checkNeighbors(bestOld, map, center);

    %% ===================================================================
    %% Future work: grid optimization
    %% ===================================================================
    % dir = best - bestOld;
    % best = optimizeGridOrigin(best, vpgon, vertices, dir);

    % Save Figure 2 frame in the animation
    writeVideo(v2, getframe(fig2));

    %%
    map = grid2D(fp0, olapx, olapy, best, roi);
    t = t + tobs;
    % t = t + tobs + slew(t, A{end-1}, A{end});
end

% End animations
close(v2)
end