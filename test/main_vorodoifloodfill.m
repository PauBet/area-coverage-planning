clc; clear all; close all;

%%
addpath(genpath(pwd));

%%
roi = [60 20;
              60 0;
              90 0;
              90 20;];

targetArea = polyshape(roi);

x = roi(:,1);
y = roi(:,2);

vertices(:,1) = x;
vertices(:,2) = y;

%%
inputkernels_galileo;
addpath('/Users/paulabetriu/Desktop/GitHub/RESSlib');
initSPICEv(fullK(METAKR));

%%
startTime = cspice_str2et('1997 DEC 16 11:15:00.000 TDB');
endTime   = cspice_str2et('1997 DEC 16 12:00:00.000 TDB');
step = 20;
inst = 'GLL_SSI';
sc = 'GALILEO ORBITER';
target = 'EUROPA';
method = 'ELLIPSOID';
theta = 0;

%%
t = startTime;
olapx = 0; olapy = 0;
[gamma(1), gamma(2)] = centroid(targetArea);
fprintc = footprint(gamma(1), gamma(2), t, inst, sc, target, theta);

% Get the footprint bounding box
bbox  = smallestBoundingBox(fprintc.bvertices(:, 1), fprintc.bvertices(:, 2));

angle = deg2rad(bbox.angle);

%gridPoints = floodFillVoronoiAlgorithm(fp.sizex, fp.sizey, olapx, olapy,...
%     gamma, areapoints, [], '8fill');
cside = closestSide(target, sc, t, roi);
tour = planSidewinderTour(target, sc, t, roi, fprintc,...
        olapx, olapy, cside, {gamma}, {});
fp = footprint(tour{1}(1), tour{2}(1), t, inst, sc, target, theta);

% Filling the region-of-interest (roi) with a footprint that is not aligned
% with the meridian-equator axes is equivalent to filling the oriented
% target area with an aligned footprint (angle = 0). Therefore, we rotate
% the region-of-interest to orient it according to the footprint
rotmat = [cos(angle)   -sin(angle);
          sin(angle)   cos(angle)];
[cx, cy] = centroid(polyshape(roi(:,1), roi(:,2)));
orientedArea  = zeros(length(roi), 2);
for j=1:length(roi)
    orientedArea(j, :) = [cx, cy]' + rotmat*(roi(j, :)' - ...
        [cx, cy]');
end
gamma = [cx, cy]' + rotmat*(gamma' - [cx, cy]');


figure
hold on; box on; grid minor;
floodFillVoronoiAlgorithm(fprintc.sizex, fprintc.sizey, t, inst, sc, target,...
    theta, olapx, olapy, tour{1}, fprintc.sizex, tour{1}, orientedArea, [], [], '8fill');
plot(polyshape(orientedArea(:,1), orientedArea(:,2)), 'FaceColor', 'none', 'EdgeColor', [0.64,0.08,0.18], ...
    'linewidth', 1.5)
xlabel('Planetocentric longitude [º]')
ylabel('Planetocentric latitude [º]')
title('ROI flood-fill discretization')
legend('Reference tile', 'Stare point')
set(gca, 'fontsize', 18)
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ROI-static-discretization','-dpdf','-r0')
hold off;

%
figure
hold on; box on; grid minor;
gamma = tour{1};
floodFillAlgorithm(fprintc.sizex, fprintc.sizey, olapx, olapy,...
     gamma, orientedArea, [], '8fill');
plot(polyshape(orientedArea(:,1), orientedArea(:,2)), 'FaceColor', 'none', 'EdgeColor', [0.64,0.08,0.18], ...
    'linewidth', 1.5)
xlabel('Planetocentric longitude [º]')
ylabel('Planetocentric latitude [º]')
title('ROI flood-fill discretization')
legend('Reference tile', 'Stare point')
set(gca, 'fontsize', 18)
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ROI-static-discretization','-dpdf','-r0')

%%
endSPICE;