clear; clc; close all;
% grid2D test

addpath(genpath(pwd))

targetArea = [ 30  20;
               30 -20;
              -30 -20;
              -30  20];

sx = 10; sy = 30;
fpref.sizex = sx;
fpref.sizey = sy;
fpref.bvertices = [ sx/2  sy/2;
                    -sx/2 sy/2;
                   -sx/2   -sy/2;
                   sx/2  -sy/2];

ang = -30;
rotmat = [cosd(ang) -sind(ang);
          sind(ang) cosd(ang)];
for i=1:length(fpref.bvertices)
    fpref.bvertices(i, :) = rotmat*fpref.bvertices(i,:)';
end

%grid2D(fpref, 0, 0, [0, 0], targetArea);

tour = planSidewinderTour('VESTA', 'DAWN', cspice_str2et('2011 SEP 30 00:00:00.00'), targetArea, fpref, ...
    [0, 0], 0, 0);

for i=1:length(tour)
    mattour(i, :) = tour{:, i};
end

plot(mattour(:,1), mattour(:,2))