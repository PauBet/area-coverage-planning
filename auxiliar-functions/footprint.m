function fp = footprint(lon, lat, instName, targetName, scName, t, method, ax, fig, video)
% [Description]

%% Pre-allocate variables
fp.bbox = []; % Minimum and maximum values of the latitude and longitude 
% bounding the footprint

% Longitude and latitude values of the footprint origin (boresight
% projection of the instrument onto the body surface)
fp.olon = lon; 
fp.olat = lat;

% Horizontal (lon) and vertical (lat) size values of the footprint
fp.sizex = 0; % longitude size [deg]
fp.sizey = 0; % latitude size [deg]

%% Parameters definition
lon = lon*cspice_rpd; % longitude from [deg] to [rad]
lat = lat*cspice_rpd; % latitude from [deg] to [rad]
[~, instFrame, boresight, bounds] = ...
    cspice_getfov(cspice_bodn2c(instName), 4); % instrument FOV's boundary
    % vectors in the instrument frame
[~, targetFrame, ~] = cspice_cnmfrm(targetName); % target frame ID in SPICE
abcorr = 'LT'; % one-way light time aberration correction parameter. 
% See cspice_spkpos for further information.
N = 200; % footprint vertices resolution

%% Instrument pointing orientation change
% Rectangular coordinates of the target point
recPoint = cspice_srfrec(cspice_bodn2c(targetName), lon, lat);
% Instrument (spacecraft) position
instpos  = cspice_spkpos(scName, t, targetFrame, abcorr, targetName);
% Distance vector to the target point from the instrument in the
% targetFrame reference system
v2 = recPoint - instpos;
% Check if the point is visible as seen from the instrument, otherwise the
% function is exited and the footprint is returned empty
if dot(v2, recPoint) > 0
    return;
end
% Rotation matrix from target to instrument frames
rotationMatrix = cspice_pxform(targetFrame, instFrame, t);
% Distance vector to the target point from the instrument in the instFrame
% reference system
v2 = rotationMatrix*v2;
% Rotation axis over which the instrument pointing has to be rotated, i.e.
% the cross vector of the boresight (origin) vector and the final pointing
% vector (v2)
rotAxis = normalize(cross(boresight, v2), 'norm');
% Phase that has to be rotated from one vector to the other
angle = cspice_vsep(boresight, v2);
% New rotation matrix that sets the new instrument pointing, i.e., towards
% the target point [lon lat]
pointingRotation = cspice_axisar(rotAxis, angle);
% New FOV bounds orientation
for i=1:length(bounds)
    bounds(:,i) = pointingRotation*bounds(:,i); % instrument FOV's boundary
    % vectors in the instrument frame (bounds)
end

%% FOV boundaries intersection with the body (footprint vertices)
% In its simplest form, the footprint should have the same number of
% vertices as boundaries has the instrument's FOV
boundPoints = zeros(3,length(bounds)); % intercept points of the FOV
% boundary, in Cartesian coordinates, and in the body-fixed reference frame
% In its simplest form, the footprint should have the same number of
% vertices as boundaries has the instrument's FOV
irr = false; % boolean to define if the instrument FOV projection is not
% enclosed in the body surface, i.e. at least one of the boundary vectors
% do not intercept the body surface
for i=1:length(bounds)
    [boundPoints(:,i), ~, ~, found] = cspice_sincpt(method, targetName, t,...
        targetFrame, abcorr, scName, instFrame, bounds(:,i));
    % If the FOV boundary does not intercept the target...
    if ~found
        irr = true;
        break;
    end
end

surfPoints = [];
count = 0;
if irr
    % For those cases where the FOV does not completely contain the target,
    % a more refined search is going to be performed in order to define the
    % limits of the footprint
    minx = min(bounds(1,:)); % minimum longitude value, in [rad]
    maxx = max(bounds(1,:)); % maximum longitude value, in [rad]
    miny = min(bounds(2,:)); % minimum latitude value, in [rad]
    maxy = max(bounds(2,:)); % maximum latitude value, in [rad]
    z = bounds(3,1); % z-coordinate of the boundary vectors
    for i=0:N
        x = (maxx - minx)*i/N + minx;
        for j=0:N
            y = (maxy - miny)*j/N + miny;
            vec = [x, y, z]';
            vec = pointingRotation*vec;
            [aux, ~, ~, found] = cspice_sincpt(method,...
                targetName, t, targetFrame, abcorr, scName,...
                instFrame, vec);
            if j == 0
                old_found = found;
                old_surfPoint = aux;
            end

            if found && (y == miny || y == maxy || x == minx || x == maxx)
                count = count + 1;
                surfPoints(:, count) = aux;
            elseif found ~= old_found
                count = count + 1;
                if old_found
                    surfPoints(:,count) = old_surfPoint;
                else
                    surfPoints(:,count)  = aux;
                end
            end
            old_found = found;
            old_surfPoint = aux;
        end
    end
else
    % The FOV projection is enclosed in the target surface so a linear
    % interpolation approximation is performed in this case (more
    % efficient)
    boundPoints(:, end+1) = boundPoints(:, 1); % close polygon
    for i=1:length(boundPoints)-1
        % linear (approximation) interpolation between vertices to define 
        % the boundary of the footprint
        v = boundPoints(:, i+1) - boundPoints(:, i); % director vector
        lambda = linspace(0, 1, N); % line parametrization
        for l=1:N
            count = count + 1;
            surfPoints(1, count) = boundPoints(1, i) + v(1)*lambda(l); % x
            surfPoints(2, count) = boundPoints(2, i) + v(2)*lambda(l); % y
            surfPoints(3, count) = boundPoints(3, i) + v(3)*lambda(l); % z
        end
    end
end

for i=1:length(surfPoints)
    [~, auxlon, auxlat] = cspice_reclat(surfPoints(:,i)); % rectangular
    % to latitudinal coordinates
    aux(i,1) = auxlon*cspice_dpr(); % longitude in [deg]
    aux(i,2) = auxlat*cspice_dpr(); % latitude in [deg]
end

% When the footprint contains the limb, the footprint intercept is
% irregular, meaning that the boundary is not a smooth curve (the limb) but
% a set of scattered points with certain deviation around the limb.
% As an approximate solution, the outer points are taken.
k = boundary(aux(:,1), aux(:,2)); % most external points

% Save footprint vertices
fp.bbox(:,1) = aux(k,1);
fp.bbox(:,2) = aux(k,2);

%% temporal
% if ~isempty(k)
%     radii = cspice_bodvrd('VESTA','RADII',3);
%     limb  = cspice_edlimb(radii(1), radii(2), radii(3), instpos);
%     a = norm(limb.semiMajor);
%     b = norm(limb.semiMinor);
%     ellipse = @(x,y) x^2/a^2 + y^2/b^2 - 1;
%     ecc = 1 - (b/a)^2;
%     theta = 0:1:360;
%     r = limb.center + limb.semiMajor*sind(theta) + limb.semiMinor*cosd(theta);
%     [~, rlon, rlat] = cspice_reclat(r);
%     
%     k = boundary(fp.bbox(:,1), fp.bbox(:,2));
%     figure
%     scatter(fp.bbox(:,1),fp.bbox(:,2))
%     hold on; box on;
%     plot(fp.bbox(k,1), fp.bbox(k,2))
%     scatter(rlon*cspice_dpr, rlat*cspice_dpr)
% 
%     figure
%     plot3(r(1,:), r(2,:), r(3,:))
% end

%% FOV representation
% This figure shows the FOV projection of the instrument over the target's
% surface, modeled as a triaxial ellipsoid. If method is equal to
% 'DSK', this option is disabled (for now, it could be considered
% to represent the DEM model in the near future, although it is going
% to be computationally demanding)

if ~isequal(method,'DSK/UNPRIORITIZED') 
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
    instpos  = cspice_spkpos(scName, t, targetFrame, 'NONE', targetName);
    rotationMatrix = cspice_pxform(instFrame, targetFrame, t);
    bbounds = zeros(3, length(bounds));
    for i=1:length(bounds)
        rot = rotationMatrix*bounds(:,i);
        bbounds(:,i) = normalize(rot, 'norm') + normalize(instpos,'norm');
    end
    pyramid_vertex = bbounds*norm(instpos);
    pyramid_vertex(:,end+1) = instpos;

    plot3(ax, instpos(1),instpos(2),instpos(3),'p','MarkerFaceColor',...
        [255 146 4]/255, 'MarkerEdgeColor', [255 146 4]/255)
    face= [2 3 5; 1 2 5; 3 4 5; 4 1 5];
    p = patch(ax, 'Faces',face,'Vertices',pyramid_vertex','Facecolor',...
        [0.66 0.85 0.41], 'FaceAlpha', 0.5);
    hold on;
    plot3(ax, surfPoints(1,:), surfPoints(2,:), surfPoints(3,:),...
        '.b','MarkerSize',4)
    view(ax, instpos)
    drawnow
    writeVideo(video, getframe(fig));
    set(p,'visible','off')
end

%% Parameters
fp.olon = lon;
fp.olat = lat;
maxlon = max(fp.bbox(:,1));
minlon = min(fp.bbox(:,1));
maxlat = max(fp.bbox(:,2));
minlat = min(fp.bbox(:,2));

mid = zeros(3,length(bounds));
for i=1:length(bounds)
    if i~=length(bounds)
        mid(:,i) = .5*(bounds(:,i) + bounds(:,i+1));
    else
        mid(:,i) = .5*(bounds(:,i) + bounds(:,1));
    end
end

fp.sizex = maxlon - minlon;
fp.sizey = maxlat - minlat;
end