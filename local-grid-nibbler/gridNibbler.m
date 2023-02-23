function A = gridNibbler(startTime, endTime, step, instName, scName, targetName, vertices, method)
%% [Description]

%% Initialization of variables
A = {};
scPos = [];
x = vertices(:,1); y = vertices(:,2);
bx = []; by = [];
poly1 = polyshape(x,y);
[~, targetFrame, ~] = cspice_cnmfrm(targetName);

%% Figure 1.
% This figure shows the FOV projection of the instrument over the target's
% surface, modeled as a triaxial ellipsoid. If method is equal to
% 'DSK', this option is disabled (for now, it could be considered
% to represent the DEM model in the near future, although it is going
% to be computationally demanding)

fig1 = figure;
ax1 = axes;
set(gcf,'units','normalized','OuterPosition',[0.0052,0.3139,0.4201,0.6532]);
hold on; grid minor; axis equal; box on; grid on;
set(gca,'Color','k','GridColor','w','MinorGridColor','w','XColor','k',...
    'YColor','k','ZColor','k','FontSize',15)
ax = gca;
ax.XAxis.TickLabelColor = 'k';
ax.YAxis.TickLabelColor = 'k';
ax.ZAxis.TickLabelColor = 'k';
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Footprint projection')
radii = cspice_bodvrd('VESTA','RADII',3);
[xe, ye, ze] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),100);
surf(xe, ye, ze, 'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.50 0.50 0.50])
quiver3(0,0,0,1.5*radii(1),0,0,'c','linewidth',2)
quiver3(0,0,0,0,1.5*radii(2),0,'g','linewidth',2)
quiver3(0,0,0,0,0,1.5*radii(3),'r','linewidth',2)

% animation figure 1
v1 = VideoWriter('footprint_projection_nibbler');
v1.FrameRate = 2;
open(v1);
writeVideo(v1,getframe(fig1));

%% Figure 2.
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex objects
ax2 = mapPlot('vesta-map.png');
fig2 = gcf;
set(gcf,'units','normalized','OuterPosition',[0.4307,0.3144,0.5656,0.6565]);
plot(ax2, polyshape(vertices), 'FaceColor', [0.93,0.69,0.13])

% animation figure 2
v2 = VideoWriter('topography_map_nibbler');
v2.FrameRate = v1.FrameRate;
open(v2);
writeVideo(v2,getframe(fig2));

%% First iteration
t = startTime;
scPos(:,end+1) = cspice_spkpos(scName, t, targetFrame, 'NONE', targetName);
% Assumption: we choose to start in one of the target's area corners,
% the one closest to the spacecraft's ground track position
[gamma(1), gamma(2)] = centroid(polyshape(vertices(:,1),vertices(:,2)));
fp0 = footprint(gamma(1), gamma(2), t, instName,  scName, targetName, 0); 
map = grid2D(fp0, 0, 0, gamma, vertices);
best = closestTargetCorner(targetName, scName, t, map);

%% Local Grid Nibbler
while ~iszero(vertices) && t <= endTime
    %
    %
    %     fpshape = polyshape(vpgon.bbox);
    %     inter = subtract(poly1, fpshape);
    %     areaI = area(inter);
    %     areaT = area(poly1);
    %     areaInter = areaT - areaI;
    %     %     if areaInter/area(fpshape) < 0.2
    %     %         newStart = closestTargetCorner(bestOld, t);
    %     %         bestNew = checkNeighbors(newStart, t);
    %     %         if score(bestNew, center) > score(best, center)
    %     %             best = newStart;
    %     %         end
    %     %     end
    
    % Footprint projection and subtraction from the ROI
    vpgon = footprint(best(1), best(2), instName, targetName, scName, ...
        t, method, [], [], []);
    poly2 = polyshape(vpgon.bbox);
    poly1 = subtract(poly1, poly2);
    vertices = poly1.Vertices;

    % Plots (figures 1 & 2)
    A{end + 1} = best;
    plot(ax2, poly2, 'FaceColor', [0.00,0.45,0.74])
    bx = [bx A{end}(:,1)];
    by = [by A{end}(:,2)];
    if length(bx) > 1
        plot(ax2, bx(end-1:end), by(end-1:end), 'w-', 'linewidth', 1)
    end
    drawnow
    if (t > startTime)
        scPos(:,end+1) = cspice_spkpos(scName, t, targetFrame, 'NONE', targetName);
        plot3(ax1, scPos(1,end-1:end),scPos(2,end-1:end),scPos(3,end-1:end),'Color', [255 146 4]/255)
    end
    groundTrack = cspice_subpnt('INTERCEPT/ELLIPSOID', targetName, t, targetFrame,...
        'NONE', scName);
    [~, sclon, sclat] = cspice_reclat(groundTrack);
    plot(ax2, sclon*cspice_dpr, sclat*cspice_dpr, 'y^', 'MarkerSize', 6)
    writeVideo(v2, getframe(fig2));

    % Starting from the previous tour point, all of its 8 neighbours are
    % going to be checked (and scored)
    bestOld = best;
    [center(1), center(2)] = centroid(polyshape(vertices(:,1), vertices(:,2)));
    best = checkNeighbors(bestOld, map, center);

    %% ===================================================================
    %% Future work: grid optimization
    %% ===================================================================
    % dir = best - bestOld;
    % best = optimizeGridOrigin(best, vpgon, vertices, dir);

    %%
    map = grid2D(fp0.sizex, fp0.sizey, fp0.overlapx, fp0.overlapy, best, vertices);
    t = t + step;
    % t = t + tobs + slew(t, A{end-1}, A{end});
end
end