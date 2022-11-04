clear all; close all; clc;


w = 5; h = 8;
targetArea = [0 0;
              0 60;
              30 100;
              60 80;
              60 0;
              30 10];

targetArea = [0 0;
              0 60;
              60 60;
              60 0];

hold on
I = flip(imread('europa-map.png')); 
imagesc([0 360],[-90 90],I);
xlim([0 360]);
ylim([-90 90]);
xlabel('longitude [°]','FontSize', 12)
ylabel('latitude [°]','FontSize', 12)

n = 10;
targetArea = add_area_manual(n);

ox = 1; oy = 2;
gridPoints = [];
gamma = [sum(targetArea(:,1)), sum(targetArea(:,2))]/length(targetArea);

figure
hold on;
xlim([min(targetArea(:,1) - 2*w) max(targetArea(:,1) + 2*w)])
ylim([min(targetArea(:,2) - 2*h) max(targetArea(:,2) + 2*h)])
plot(polyshape(targetArea))
tic
floodFillAlgorithm(w, h, ox, oy, gamma, targetArea, gridPoints, '4fill')
%[gridx, gridy] = grid2D(w, h, ox, oy, gamma, targetArea);
toc