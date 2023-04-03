clear all; close all; clc;

gamma = [0 0];
w = 10; h = 6;
fpx = [gamma(1) - w/2, gamma(1) - w/2, gamma(1) + w/2, gamma(1) + w/2];
fpy = [gamma(2) + h/2, gamma(2) - h/2, gamma(2) - h/2, gamma(2) + h/2];
angle = 0:15:360;

rotfp = zeros(length(fpx), 2);
for i=1:length(angle)
    rotmat = [cosd(angle(i))   -sind(angle(i));
        sind(angle(i))   cosd(angle(i))];
    for j=1:length(fpx)
        rotfp(j, :) = gamma' + rotmat*([fpx(j), fpy(j)] - gamma)';
    end
    figure
    hold on; box on; grid minor;
    plot(polyshape(rotfp(:, 1), rotfp(:, 2)));
    set(gca, 'fontsize', 15)
    axis equal

    bbox = smallestBoundingBox(rotfp(:, 1), rotfp(:, 2));
    disp("input angle: " + angle(i) + " ; output angle: " + bbox.angle);
end