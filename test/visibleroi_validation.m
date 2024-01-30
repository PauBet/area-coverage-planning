clc; close all; clear all;
% validation of visibleroi function

% SPICE initialization
addpath(genpath(pwd));
resslib = '/Users/paulabetriu/Desktop/GitHub/RESSlib'; % library of
% SPICE functions (specially relevant for SPICE initialization in this
% program)
addpath(resslib);
run(fullfile(pwd, 'input', 'galileo', 'inputkernels.m'));
initSPICEv(fullK(METAKR));

%% ROI 1: Does not intercept the a.m. (but the limb does)
% Definition of parameters
roi = [-55   20;
       -85   20;
       -85  -20;
       -55  -20;];

target = 'EUROPA';
obs    = 'GALILEO ORBITER';
[~, targetframe, ~] = cspice_cnmfrm(target); % body-fixed frame
% time window
inittime = cspice_str2et('1998 MAR 29 13:30:00.000 TDB');
stoptime = cspice_str2et('1998 MAR 29 15:30:00.000 TDB');
et = inittime:10*60:stoptime;

% call function
for i=1:length(et)
    % compute ground track
    sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
        et(i), targetframe, 'NONE', obs);
    [~, sclon, sclat] = cspice_reclat(sctrack);
    figure
    plot(polyshape(roi(:, 1), roi(:, 2)), 'FaceColor', 'none', 'edgecolor', 'r')
    hold on; grid minor; box on; axis tight;
    [vroi, poly1] = visibleroi(roi, et(i), target, obs);
    plot(polyshape(vroi(:, 1), vroi(:, 2)))
    plot(poly1, 'facecolor', 'none', 'edgecolor', 'b')
    plot(sclon*cspice_dpr, sclat*cspice_dpr, 'b^', 'markersize', 14)
    xlim([-180 180])
    ylim([-90 90])
end

%% ROI 2: Both the ROI and the limb intercept the anti-meridian line
roi = [-177  3;
       -177 -3;
        177 -3;
        177  3];

% call function
[sroi(:, 1), sroi(:, 2)] = amsplit(roi(:, 1), roi(:, 2));
for i=1:length(et)
    % compute ground track
    sctrack = cspice_subpnt('INTERCEPT/ELLIPSOID', target,...
        et(i), targetframe, 'NONE', obs);
    [~, sclon, sclat] = cspice_reclat(sctrack);
    figure
    plot(polyshape(sroi(:, 1), sroi(:, 2)), 'FaceColor', 'none', 'edgecolor', 'r')
    hold on; grid minor; box on; axis tight;
    [vroi, poly1] = visibleroi(roi, et(i), target, obs);
    plot(polyshape(vroi(:, 1), vroi(:, 2)))
    plot(poly1, 'facecolor', 'none', 'edgecolor', 'b')
    plot(sclon*cspice_dpr, sclat*cspice_dpr, 'b^', 'markersize', 14)
    xlim([-180 180])
    ylim([-90 90])
end