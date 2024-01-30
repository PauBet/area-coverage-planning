clc; clear; close all;
% Generate random points
% n = 50;
% points = rand(2, n);
% 
% x = linspace(0, 1, 20);
% y = linspace(0, 1, 20);
% x = randi([0, 20], [1, 20]);
% y = randi([0, 20], [1, 20]);
% x = x/20; y = y/20;
% X = [];
% for i=1:length(x)
%     for j=1:length(y)
%     point = [x(i) y(j)];
%     X = [X; point];
%     end
% end
% 
% % Compute Voronoi diagram
% %voronoiDiagram = voronoi(x, y);
% [vv, vc] = voronoin(X);
% 
% % % Extract the finite Voronoi edges
% % finiteEdges = voronoiDiagram(1, :);
% % for i = 2:size(voronoiDiagram, 1)
% %     if ~any(voronoiDiagram(i, :) == Inf)
% %         finiteEdges = [finiteEdges; voronoiDiagram(i, :)];
% %     end
% % end
% 
% [x, y] = meshgrid(x, y);
% 
% % Plot Voronoi diagram
% figure
% hold on
% plot(x, y, 'ko');
% vv(vv == inf) = 1e10;
% for i=1:length(vc)
%     ind = vc{i};
%     plot(polyshape(vv(ind, 1), vv(ind, 2)))
% end

% 1-D
n = 30;
xmax = 40;
x = linspace(0, xmax, n);
y = zeros(1, n);
refx = xmax/n; refy = 1;
sizex = refx*randi([9, 11], [1, length(x)]);
sizex = sizex./10;
w = refx; h = refy;
sizex(1) = refx;

figure
hold on;
corrx = linspace(0, xmax, n);
gamma = [corrx(1) 0];
pp = polyshape([gamma(1)-w/2, gamma(1)-w/2, gamma(1)+w/2, gamma(1) + w/2, gamma(1)-w/2],[gamma(2)+ h/2, gamma(2)- h/2, gamma(2)- h/2, gamma(2)+ h/2, gamma(2)+ h/2]);
plot(pp, 'FaceColor', [0.93,0.69,0.13], 'FaceAlpha', 0.2);
for i=2:length(corrx)
    s1 = sizex(i-1);
    s2 = sizex(i);
    deltax = corrx(i) - corrx(i-1) - .5*(s1 + s2);
    corrx(i) = corrx(i) - deltax;
    w = s2;
    gamma = [corrx(i) 0];
    pp = polyshape([gamma(1)-w/2, gamma(1)-w/2, gamma(1)+w/2, gamma(1) + w/2, gamma(1)-w/2],[gamma(2)+ h/2, gamma(2)- h/2, gamma(2)- h/2, gamma(2)+ h/2, gamma(2)+ h/2]);
    plot(pp, 'FaceColor', [0.93,0.69,0.13], 'FaceAlpha', 0.2);
end
plot(x, y, 'r^')
hold on; grid minor; axis tight;
plot(corrx, y, 'bp')
set(gca, 'fontsize', 20)

%
n = 10;
x = linspace(0, n, n);
y = linspace(0, n, n);
[x, y] = meshgrid(x, y);
refx = 1; refy = 1;

%
sizex = randi([5, 15], [1, length(x)]);
sizey = randi([5, 15], [1, length(x)]);
sizex = sizex./15;
sizey = sizey./15;

corrx = linspace(0, n, n);
corry = linspace(0, n, n);
auxx = corrx;
auxy = corry;
for i=2:length(corrx)
    scalex = sizex(i) / refx;
    spacex = refx * scalex;
    corrx(i)  = auxx(i-1) + spacex;
end
for j=2:length(corry)
    scaley = sizey(j) / refy;
    spacey = refy * scaley;
    corry(j)  = auxy(j-1) + spacey;
end
[x1, y1] = meshgrid(corrx, corry);

figure
hold on;
for i=2:size(x1,1)
    for j=2:size(x1,2)
        w = x1(i, j) - x1(i, j-1);
        h = y1(i, j) - y1(i-1, j);
        gamma = [x1(i, j) y1(i, j)];
        pp = polyshape([gamma(1)-w/2, gamma(1)-w/2, gamma(1)+w/2, gamma(1) + w/2, gamma(1)-w/2],[gamma(2)+ h/2, gamma(2)- h/2, gamma(2)- h/2, gamma(2)+ h/2, gamma(2)+ h/2]);
        plot(pp, 'FaceColor', [0.93,0.69,0.13], 'FaceAlpha', 0.2);
        drawnow
    end
end
plot(x, y, 'r^')
hold on; grid minor; axis tight;
plot(x1, y1, 'bp')
set(gca, 'fontsize', 20)


axis equal