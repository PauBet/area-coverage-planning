function A = replanningSidewinder(startTime, endTime, tobs, instName, ...
    scName, targetName, vertices, olapx, olapy, method)
% Optimization of the Sidewinder algorithm, adapted from [1]. 
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        A = replanningSidewinder(startTime, endTime, tobs, instName, ...
%               scName, targetName, vertices, method)
%
% Inputs:
%   > startTime:    start time of the planning horizon, in TDB seconds past
%                   J2000 epoch
%   > endTime:      end time of the planning horizon, in TBD seconds past
%                   J200 epoch
%   > tobs:         observation time, i.e. the minimum time that the 
%                   instrument needs to perform an observation, in seconds
%   > instName:     string name of the instrument
%   > scName:       string name of the spacecraft
%   > targetName:   string name of the target body
%   > vertices:     matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D. 
%       # vertices(:,1) correspond to the x values of the vertices
%       # vertices(:,2) correspond to the y values of the vertices
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in deg
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in deg
%   > method:       string defining the computation method needed for the
%                   cspice_sincpt function. See cspice_sincpt function
%                   description for more information
%       Possible values:
%       * 'ELLIPSOID' builds an ellipsoid to model the body surface
%       * 'DSK/UNPRIORITIZED[/SURFACES = <surface list>]' uses topographic
%       data to model the body surface
% 
% Outputs:
%   > A:            cell matrix of the successive instrument observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

%%
% Pre-allocate variables
A = {}; % List of observations (successive boresight ground track position)
tour = {}; % Sorted ROI's grid discretization (Boustrophedon decomposition)

% Define target area as a polygon
x = vertices(:,1); y = vertices(:,2);
poly1 = polyshape(x,y);

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
radii = cspice_bodvrd(targetName,'RADII',3);
[xe, ye, ze] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),100);
surf(xe, ye, ze, 'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.50 0.50 0.50])
quiver3(0,0,0,1.5*radii(1),0,0,'c','linewidth',2)
quiver3(0,0,0,0,1.5*radii(2),0,'g','linewidth',2)
quiver3(0,0,0,0,0,1.5*radii(3),'r','linewidth',2)

% animation figure 1
v1 = VideoWriter('footprint_projection_rsidewinder');
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
v2 = VideoWriter('topography_map_rsidewinder');
v2.FrameRate = v1.FrameRate;
open(v2);
writeVideo(v2,getframe(fig2));

%% Replanning Sidewinder algorithm
% The first time iteration is the starting time in the planning horizon
t = startTime;

% Compute spacecraft position in the body-fixed frame
[~, targetFrame, ~] = cspice_cnmfrm(targetName); % body-fixed frame
scPos = cspice_spkpos(scName, t, targetFrame, 'NONE', targetName);

% Closest polygon side to the spacecraft's ground track position
closestSide = getClosestSide(targetName, scName, t, vertices);

while ~iszero(vertices) && t <= endTime

    if t == startTime
        % Initial 2D grid layout discretization: the instrument's FOV is
        % going to be projected onto the uncovered area's centroid and the
        % resulting footprint shape is used to set the grid spatial
        % resolution
        [gamma(1), gamma(2)] = centroid(polyshape(vertices(:,1),...
            vertices(:,2)));
        fprint0 = footprint(gamma(1), gamma(2), instName, targetName,...
            scName, t, method, ax1, [], v1);  % centroid footprint
    else
        % For the subsequent time steps, the following observation
        % programmed in the previous tour calculation is used to discretize
        % the remaining uncovered area
        gamma = tour{1};
        fprint0 = footprint(gamma(1), gamma(2), instName, targetName, ...
            scName, t, method, ax1, [], v1);
    end
    
    % Grid optimization
    if isempty(tour)
        [gamma(1), gamma(2)] = centroid(polyshape(vertices(:,1),...
            vertices(:,2)));
    else
        % Direction of the next observation in the coverage path
        dir = gamma - A{end};
        % Optimization of the grid in order to avoid gridlock
        gamma = optimizeGridOrigin(gamma, fprint0, olapx, olapy,...
            vertices, dir);
    end

    % Sorted list of grid points according to the sweeping/coverage path
    % (see Boustrophedon decomposition)
    tour = planSidewinderTour(closestSide, vertices, fprint0, gamma,...
        olapx, olapy);

    % Compute the footprint of each point in the tour successively and
    % subtract the corresponding area from the target polygon
    a = tour{1}; % observation
    tour(1) = []; % delete this observation from the planned tour
    A{end + 1} = a; % add it in the list of planned observations
    fprinti = footprint(a(1), a(2), instName, targetName, scName, ...
        t, method, ax1, fig1, v1); % compute the observation's footprint
    if ~isempty(fprinti.bbox)
        poly2 = polyshape(fprinti.bbox); % create footprint polygon
        poly1 = subtract(poly1, poly2); % update uncovered area
        vertices = poly1.Vertices;

        %% Plots

        % Footprint plot in Figure 2
        plot(ax2, poly2, 'FaceColor', [0.00,0.45,0.74])
        if length(A) > 1
            plot(ax2, [A{end-1}(1) A{end}(1)], [A{end-1}(2) A{end}(2)],...
                'w-', 'linewidth', 1)
        end
        drawnow
        
        % Spacecraft position plot in Figure 1
        if (t > startTime)
            scPos(:,end+1)  = cspice_spkpos(scName, t, targetFrame,...
                'NONE', targetName);
            plot3(ax1, scPos(1,end-1:end),scPos(2,end-1:end),...
                scPos(3,end-1:end),'Color', [255 146 4]/255)
        end

        % Spacecraft ground track position in Figure 2
        groundTrack = cspice_subpnt('INTERCEPT/ELLIPSOID', targetName,...
            t, targetFrame,...
            'NONE', scName);
        [~, sclon, sclat] = cspice_reclat(groundTrack);
        plot(sclon*cspice_dpr, sclat*cspice_dpr, 'y^', 'MarkerSize', 6)

        % Save Figure 2 frame in the animation
        writeVideo(v2, getframe(fig2));

        % New time iteration
        t = t + tobs;
        %t = t + tobs + slewDur(t, A{end}, A{end + 1}); % future work
    end
end
% End animations
close(v1)
close(v2)
end