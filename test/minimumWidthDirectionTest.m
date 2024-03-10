% Test minimumWidthDirection
clc; close all; clear all;

roi = [10  5;
       10  0;
       20  0;
       20  5];

poly1 = polyshape(roi(:, 1), roi(:, 2));
[cx, cy] = centroid(poly1);

figure
hold on; box on; grid minor; axis equal;
plot(poly1);

fprintf('Minimum Width Direction (thetamin): %.2fº\n', minimumWidthDirection(roi(:, 1), roi(:, 2)));

rotmat = @(angle) [cosd(angle) -sind(angle);
                   sind(angle)  cosd(angle)];

for i=1:size(roi, 1)
    rroi(i, :) = rotmat(45)*(roi(i, :)' - [cx, cy]');
end
poly2 = polyshape(rroi(:, 1), rroi(:, 2));

figure
hold on; box on; grid minor; axis equal;
plot(poly2);

fprintf('Minimum Width Direction (thetamin): %.2fº\n', minimumWidthDirection(rroi(:, 1), rroi(:, 2)));

figure
hold on; grid minor; box on;
xlim([-180 180])
ylim([-90 90])
n = 10;
for i=1:n
    roi(i, :) = ginput(1);
    plot(roi(i, 1) , roi(i, 2),'.','MarkerSize',8,'Color','blue')
    hold on;
end
[roi(:, 1), roi(:, 2)] = sortcw(roi(:, 1), roi(:, 2));
poly3 = polyshape(roi(:, 1), roi(:, 2));
plot(poly3)
[cx, cy] = centroid(poly3);
fprintf('Minimum Width Direction (thetamin): %.2fº\n', minimumWidthDirection(roi(:, 1), roi(:, 2)));