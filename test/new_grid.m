clc; close all; clear all;
% Revision of Sidewinder:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see

% Load mission info (kernels, SPICE ids, etc.)
input_data;
inittime = cspice_str2et('1998 MAR 29 12:53:00.000 TDB'); % closest approach
tcadence = 10; % observation cadence [s]

% Define ROI (Annwn Regio)
roi = [55 30;
       55 15;
       55 10;
       35 10;
       35 30;];

% Point camera at ROI's centroid
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));

% Get camera's FOV boundaries and boresight (when pointing at centroid)
[fovbounds, boresight, rotmat] = instpointing(inst, target, sc, inittime, cx, cy);

% Build focal plane
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
vertex = cspice_spkpos(sc, inittime, targetframe, 'NONE', target);
point = vertex + fovbounds(:, 1);
plane = cspice_nvp2pl(boresight, point);

% Intersect ROI with focal plane
body = cspice_bodn2c(target);
spoint = zeros(size(roi, 1), 3);
for i=1:size(roi, 1)
    dir = -trgobsvec(roi(i, :), inittime, target, sc);
    [found, spoint(i, :)] = cspice_inrypl(vertex, dir, plane);
    if found == 0
        disp("No intersection");
    end
end
dir = normalize(dir, 'norm');
% figure
% plot3(spoint(:, 1), spoint(:, 2), spoint(:, 3), 'g*')
% hold on; box on;
% quiver3(vertex(1), vertex(2), vertex(3), boresight(1), boresight(2), boresight(3))
% quiver3(vertex(1), vertex(2), vertex(3), dir(1), dir(2), dir(3))
% points = vertex + fovbounds;
% fill3(points(1, :), points(2, :), points(3, :), 'r')

% Build grid 2D in the focal plane
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
[~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
fpref.width = width;
fpref.height = height;
fpref.angle = 0;
ovlapx = 20; ovlapy = 20;
tArea = zeros(length(spoint), 3);
for i=1:length(spoint)
    vpoint = -(vertex - spoint(i, :)');
    tArea(i, :) = rotmat\vpoint;
end
targetArea = tArea(:, 1:2);

% Focal plane grid discretization
grid = grid2D(fpref, ovlapx, ovlapy, [0, 0], targetArea);

% figure
plot(polyshape(tArea(:, 1), tArea(:, 2)))
hold on;
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        sp = grid{i, j};
        if ~isempty(sp)
            plot(sp(1), sp(2), 'b^')
        end
    end
end

% % Intersect each tile with the target's surface
% method = 'ELLIPSOID'; % assumption: ray intercept function is going to
% % model the target body as a tri-axial ellipsoid
% count = 0;
% fplist = footprint(inittime, inst, sc, target, 'lowres', cx, cy);
% for i=1:size(grid, 1)
%     for j=1:size(grid, 2)
%         sp = grid{i, j};
%         if ~isempty(sp)
%             % p = zeros(4, 3);
%             % 
%             % % Build tile
%             % p(1, 1:2) = sp + [ width/2,   height/2];
%             % p(2, 1:2) = sp + [-width/2,  height/2];
%             % p(3, 1:2) = sp + [-width/2, -height/2];
%             % p(4, 1:2) = sp + [ width/2,  -height/2];
%             % tile = polyshape(p(:, 1), p(:, 2));
%             % 
%             % % 3D coordinates
%             % p(:, 3) = 1;
%             % 
%             % % Transform coordinates to targetframe
%             % p_body = zeros(4, 3);
%             % for k=1:4
%             %     aux = rotmat*p(k, :);
%             %     p_body(k, :) = vertex - aux;
%             % end
%             % 
%             % % Intercept tile with body surface
%             % flag = false;
%             % for k=1:4
%             %     [xpoint(k, :), ~, ~, found] = cspice_sincpt(method, target, t,...
%             %         targetframe, 'NONE', sc, p_body(k, :), targetframe);
%             %     if ~found
%             %         flag = true;
%             %         break;
%             %     end
%             % end
%             % if flag
%             %     gridt{i, j} = [];
%             % else
%             %     fp = footprint(t, inst, sc, target, 'lowres', lon, lat);
%             %     gridt{i, j} = xpoint;
%             % end
% 
%             p = zeros(3, 1);
%             p(1:2) = sp;
%             p(3) = 1;
%             p_body = rotmat*p;
%             [xpoint, ~, ~, found] = cspice_sincpt(method, target, inittime,...
%                 targetframe, 'NONE', sc, targetframe, p_body);
%             if found
%                 count = count + 1;
%                 [~, lon, lat] = cspice_reclat(xpoint);
%                 fplist(count) = footprint(inittime, inst, sc, target, 'lowres', lon*cspice_dpr, lat*cspice_dpr);
%             else
%             end
%         end
%     end
% end

% Boustrophedon decomposition
tour = boustrophedon(grid, 'south', 'east');

% Intersect each tile with the target's surface
method = 'ELLIPSOID'; % assumption: ray intercept function is going to
% model the target body as a tri-axial ellipsoid
count = 0;
fplist = footprint(inittime, inst, sc, target, 'lowres', cx, cy);
t = inittime;
for i=1:length(tour)
    sp = tour{i};
    if ~isempty(sp)
        p = zeros(3, 1);
        p(1:2) = sp;
        p(3) = 1;
        p_body = rotmat*p;
        [xpoint, ~, ~, found] = cspice_sincpt(method, target, inittime,...
            targetframe, 'NONE', sc, targetframe, p_body);
        if found
            count = count + 1;
            [~, lon, lat] = cspice_reclat(xpoint);
            fplist(count) = footprint(inittime, inst, sc, target, 'lowres', lon*cspice_dpr, lat*cspice_dpr);
        else
        end
        t = t + tcadence;
    end
end

% Plot footprints
figure
plot(polyshape(roi(:, 1), roi(:, 2)))
box on;
for i=1:length(fplist)
    hold on;
    plot(polyshape(fplist(i).bvertices(:, 1), fplist(i).bvertices(:, 2)))
    plot(fplist(i).olon, fplist(i).olat, 'b^')
end
   
% Test function with limb
inittime = cspice_str2et('1998 MAR 29 14:02:00.000 TDB'); % closest approach
roi = [-55  20;
       -95  20;
       -95 -20;
       -55 -20];
[roi, polyroi] = visibleroi(roi, inittime, target, sc);
newroi = interppolygon(roi);
tour = planSidewinderTour2(target, newroi, sc, inst, inittime, ovlapx, ovlapy);

figure
plot(polyroi)
box on;
t = inittime;
for i=1:length(tour)
    hold on;
    lon = tour{i}(1); lat = tour{i}(2);
    fp = footprint(t, inst, sc, target, 'lowres', lon, lat);
    plot(polyshape(fp.bvertices(:, 1), fp.bvertices(:, 2)))
    t = t + tcadence;
end