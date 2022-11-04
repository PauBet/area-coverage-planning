clc; close all; clear all;
% Paula Betriu - August 2022
% Approximation Planning Algorithms
% This program implements a series of area coverage algorithms, based on 
% [1]. These algorithms aim to built the optimal plan for an instrument
% to represent or cover a specific area from a body surface.
%
% [1] Shao E., Byon A., Davies C., Davis E., Knight R., Lewellen G., 
% Trowbridge M. and Chien S. Area Coverage Planning with 3-axis Steerable,
% 2D Framing Sensors. 2018.

%%
addpath(genpath(pwd));

%%
%mapPlot('vesta-map.png');

% areapoints = [10 -10;
%               10 -50;
%               60 -60;
%               60 -10];
% 
% areapoints = [10 -15;
%               60 -15;
%               60 0 ;
%               10 0];

areapoints = [10 -10;
              10 -50;
              60 -60;
              60 -10;];

targetArea = polyshape(areapoints);

x = areapoints(:,1);
y = areapoints(:,2);

vertices(:,1) = x;
vertices(:,2) = y;

%%
inputkernels;
addpath('C:\Users\roure\Desktop\science-optimizer\RESSlib\');
initSPICEv(fullK(METAKR));
cspice_furnsh('C:\Users\roure\Desktop\science-optimizer\kernels\dsk\vesta512.bds')

%%
startTime = cspice_str2et('2011 SEP 30 2:00:00.000 TDB');
endTime   = cspice_str2et('2011 SEP 30 2:30:00.000 TDB');
step = 60;
instName = 'DAWN_FC2';
scName = 'DAWN';
targetName = 'VESTA';
method = 'ELLIPSOID';
A = sidewinder(startTime, endTime, step, instName, scName, targetName,...
vertices, 2, 2, method);
%A = replanningSidewinder(startTime, endTime, step, instName, scName, targetName, vertices, 2, 2, method);
%A = gridNibbler(startTime, endTime, step, instName, scName, targetName, vertices, method);

%%
endSPICE;