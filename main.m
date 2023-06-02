clc; clear all; close all
% Paula Betriu - August 2022
% Approximation Planning Algorithms
% This program implements a series of area coverage algorithms, based on 
% [1]. These algorithms aim to build the optimal plan for an instrument
% to represent or cover a specific area from a body surface.
%
% [1] Shao E., Byon A., Davies C., Davis E., Knight R., Lewellen G., 
% Trowbridge M. and Chien S. Area Coverage Planning with 3-axis Steerable,
% 2D Framing Sensors. 2018.

% Load mission info (kernels, SPICE ids, etc.)
input_data;

% Choose mosaic algorithm: 'sidewinder', 'r_sidewinder', 'onlinefrontier',
% 'gridnibbler'
tilealg = 'r_sidewinder';

% File log...
fid = fopen(fullfile(pwd, [tilealg, '.txt']), 'w');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid, '%s \n\nMosaic algorithm: %s\n', datetime("now"), tilealg);
fprintf(fid, 'Mission: %s\n', sc);

% Pre-allocation of variables... 
stoptime = cspice_str2et('1998 MAY 30 00:00:00.000 TDB'); % mosaic end (max)
tobs  = 10; % [s] between observations
olapx = 10; % [%] of overlap in x direction
olapy = 10; % [%] of overlap in y direction
videosave = 1; % =1 saves a video with the evolution of the coverage map
% along the observation plan

fprintf(fid, 'Obs time [s]: %.1f\n', tobs);
fprintf(fid, 'Overlap x [%%]: %.1f\n', olapx);
fprintf(fid, 'Overlap y [%%]: %.1f\n', olapy);

% Regions of interest
% Annwn Regio [lon, lat] = [40, 20]º
roi{1} = [50 30;
          50 15;
          30 15;
          30 30;];

inittime{1} = cspice_str2et('1998 MAR 29 12:58:00.000 TDB'); % closest approach
regions{1} = 'Annwn Regio';
    
% Pwyll crater [lon, lat] = [90, -25]º;
roi{2} = [ 78   -18;
           90   -14;
          109   -13;
          101   -23;
           99   -30;
           85   -33;
           71   -28;];

inittime{2} = cspice_str2et('1998 MAR 29 12:38:00.000 TDB'); % closest approach
regions{2} = 'Pwyll Crater';

% Tara Regio [lon, lat] = [-75, -10]º;
roi{3} = [-55   20;
          -85   20;
          -85  -20;
          -55  -20;]; % roi of roi polygon

inittime{3} = cspice_str2et('1998 MAR 29 14:04:00.000 TDB'); % closest approach
regions{3} = 'Tara Regio';

% Cilix crater [lon, lat] = [180, 0]º;
roi{4} = [-177  3;
          -177 -3;
           177 -3;
           177  3;];

inittime{4} = cspice_str2et('1998 MAR 29 13:40:00.000 TDB'); % closest approach
regions{4} = 'Cilix Crater';

% Taliesin [lon, lat] = [-138, -23]º;
roi{5} = [-160  -10;
          -160  -30;
          -130  -30;
          -130  -10;]; % roi of roi polygon

inittime{5} = cspice_str2et('1998 MAR 29 14:21:00.000 TDB');
regions{5} = 'Taliesin';

% Niamh
roi{6} = [150  25;
          150  15;
          135  15;
          135  25;]; % roi of roi polygon

inittime{6} = cspice_str2et('1998 MAR 29 13:29:00.000 TDB'); % closest approach
regions{6} = 'Niamh';

% Sort by time
aux = cell2mat(inittime);
[~, idxtime] = sort(aux);
regions  = regions(idxtime);
inittime = inittime(idxtime);
roi = roi(idxtime);

% Save ROI in log file
fprintf(fid, '\n============ INPUT ============\n');
for i=1:length(regions)
    fprintf(fid, '\n************ ROI %d ************\n', i);
    fprintf(fid, 'Name: %s\n', regions{i});
    fprintf(fid, 'Polygon:\n');
    fprintf(fid, '\t%s     \t%s\n', 'lon', 'lat');
    for j=1:size(roi{i}, 1)
        fprintf(fid, '\t%+d     \t%+d\n', roi{i}(j, 1), roi{i}(j, 2));
    end
    fprintf(fid, 'Observation time: %s\n', cspice_et2utc(inittime{i}, 'C', 0));
end

% Compute mosaics
fprintf(fid, '\n============ OUTPUT ============\n');

% Coverage figure:
% This figure plots the FOV footprint in a 2D topography map of the target 
% body. This can only be enabled for convex planetary bodies
ax = mapPlot('europa-map-shifted.png');
c  = {'g', 'y', 'c', 'm', 'r', 'k'};
% Define the RdYlBu colormap
map = getPyPlot_cMap('tab20', 7, [], '/usr/local/bin/python3');
for i=1:length(regions)
    [x, y] = amsplit(roi{i}(:, 1), roi{i}(:, 2));
    plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
        map(i, :), 'linewidth', 2)
end
legend(regions, 'AutoUpdate', 'off')

%% Mosaic algorithms
%v2 = VideoWriter('topography_map_sidewinder', 'MPEG-4');
%v2.FrameRate = 2;
%open(v2)
%writeVideo(v2, getframe(gcf));
obsDic = dictionary();
for i=1:1
    %     if abs(min(roi{i}(:, 1)) -  max(roi{i}(:, 1))) > 180
%         xlim([175  180])
%         ylim([min(roi{i}(:, 2)) - 10  max(roi{i}(:, 2)) + 10])
%     else
%         xlim([min(roi{i}(:, 1)) - 10  max(roi{i}(:, 1)) + 10])
%         ylim([min(roi{i}(:, 2)) - 10  max(roi{i}(:, 2)) + 10])
%     end

    switch tilealg
        case 'sidewinder'
            [A, cv, fplist] = sidewinder(inittime{i}, stoptime, tobs, ...
                inst, sc, target, roi{i}, olapx, olapy, ax, map(i, :), []);

        case 'r_sidewinder'
            [A, cv, fplist] = replanningSidewinder(inittime{i}, ...
                stoptime, tobs, inst, sc, target, roi{i}, olapx, olapy, ...
                ax, []);

        case 'onlinefrontier'
            [A, cv, fplist] = frontierRepair(inittime{i}, stoptime, ...
                tobs, inst, sc, target, roi{i}, olapx, olapy, ax, c{i}, 0);

        case 'gridnibbler'
            A = gridNibbler(startTime, endTime, step, instName, scName, ...
                targetName, roi, 20, 20, 1);
    end

    clear grid2D planSidewinderTour;
    obsDic(regions{i}) = {fplist};

    % Save summary of results in log file
    fprintf(fid, '\n************ ROI %d ************\n', i);
    fprintf(fid, 'Coverage [%%]: %.3f\n', cv*100);
    fprintf(fid, 'Make-span [s]: %.3f \n', fplist(end).t - fplist(1).t);
    fprintf(fid, 'Nº footprint: %d \n', (fplist(end).t - fplist(1).t)/tobs + 1);

    % Re-plot the ROI (for aesthetic purposes)
    [x, y] = amsplit(roi{i}(:, 1), roi{i}(:, 2));
    plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
        map(i, :), 'linewidth', 1)
end
set(gca, 'xlim', [-180 180], 'ylim', [-90 90])
%writeVideo(v2, getframe(gcf));
%close(v2);

% Analyze metrics
emnang    = zeros(1, length(fplist));
illzntang = zeros(1, length(fplist));
phsang    = zeros(1, length(fplist));
et        = zeros(1, length(fplist));
obsvec    = zeros(3, length(fplist));
illvec    = zeros(3, length(fplist));
meanres   = zeros(1, length(fplist));
ifov = 10e-6;  % ifov resolution [rad/px], retrieved from [Carr1995]
for j=1:length(regions)
    fplist = obsDic(regions{j});
    fplist = fplist{1};
    for i=1:length(fplist)
        et(i) = fplist(i).t;
        srfpoint = [fplist(i).olon, fplist(i).olat];
        emnang(i)    = emissionang(srfpoint, et(i), fplist(i).target, ...
            fplist(i).sc);
        illzntang(i) = illzenithang(srfpoint, et(i), fplist(i).target);
        phsang(i)    = phaseang(srfpoint, et(i), fplist(i).target, ...
            fplist(i).sc);
        obsvec(:, i) = trgobsvec(srfpoint, et(i), fplist(i).target, fplist(i).sc);
        illvec(:, i) = trgillvec(srfpoint, et(i), fplist(i).target);

        % compute resolution
        meanres(i) = pointres(ifov, srfpoint, et(i), fplist(i).target, ...
            fplist(i).sc);

        fplist(i).res = meanres(i);
    end
    obsDic(regions{j}) = {fplist};
end

figure
plot(et, emnang, '+:', 'linewidth', 1)
hold on; box on; axis tight; grid minor;
plot(et, illzntang, '+:', 'linewidth', 1)
plot(et, phsang, '+:', 'linewidth', 1)
legend('emission angle', 'illumination angle', 'phase angle', 'location', 'best')
xlabel('Time [seconds past J2000]')
ylabel('Angle [deg]')
set(gca, 'fontsize', 18)

figure
plot(et, vecnorm(obsvec), '+:', 'linewidth', 1)
hold on; box on; axis tight; grid minor;
xlabel('Time [seconds past J2000]')
ylabel('Distance [km]')
title('Surface point - observer distance')
set(gca, 'fontsize', 18)

figure
plot(et, vecnorm(illvec), '+:', 'linewidth', 1)
hold on; box on; axis tight; grid minor;
xlabel('Time [seconds past J2000]')
ylabel('Distance [km]')
title('Surface point - Sun distance')
set(gca, 'fontsize', 18)

figure
plot(et, meanres, '+:', 'linewidth', 1)
hold on; box on; axis tight; grid minor;
title('Resolution evolution')
xlabel('Time [seconds past J2000]')
ylabel('Resolution [km/px]')
set(gca, 'fontsize', 18)

% Footprint resolution visualization
maxres = -inf;
minres = inf;
res = cell(1, 4);
for j=1:length(regions)
    A = obsDic(regions{j});
    fplist = A{1};
    res{j} = zeros(1, length(fplist));
    for i=1:length(fplist)
        res{j}(i) = fplist(i).res;
        if res{j}(i) > maxres
            maxres = res{j}(i);
        end
        if res{j}(i) < minres
            minres = res{j}(i);
        end
    end
end
maxres = maxres*1.001;
N = 20;
map = getPyPlot_cMap('RdYlBu_r', N, [], '/usr/local/bin/python3');
c = map;
resrange = maxres - minres;
delta = resrange / N;
ax = mapPlot('europa-map-shifted.png');
grid minor; box on;
for rr=1:length(regions)
    A = obsDic(regions{rr});
    fplist = A{1};
    s0 = minres;
    for i=1:N
        s = i*delta + minres;
        cond = res{rr} >= s0 & res{rr} < s;
        if i==N
            cond = res{rr} >= s0 & res{rr} <= s;
        end
        ind = find(cond);
        for k=1:length(ind)
            j = ind(k);
            hold on;
            poly2 = polyshape(fplist(j).bvertices(:, 1), ...
                fplist(j).bvertices(:, 2));
            cc = c(i, :);
            plot(ax, poly2, 'FaceColor', cc, 'EdgeColor', ...
                cc, 'linewidth', 1, 'FaceAlpha', 0.2)
        end
        s0 = s;
    end
end
colormap(map)
c = colorbar();
c.Label.String = 'Resolution [km/px]';
c.Ticks = linspace(minres, maxres, 10);
clim([minres maxres])
title('Coverage map and resolution')
for i=1:length(regions)
    [x, y] = amsplit(roi{i}(:, 1), roi{i}(:, 2));
    plot(ax, polyshape(x, y), 'FaceColor', 'none', 'EdgeColor', ...
        'w', 'linewidth', 1.5)
end
set(gca, 'fontsize', 18)

% End log file
fclose(fid);

% End SPICE
endSPICE;