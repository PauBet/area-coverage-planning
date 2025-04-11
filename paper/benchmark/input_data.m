%% Relevant paths
kernelpath = '/Users/paulabetriu/Desktop/GitHub/area-coverage-planning/input';
resslib    = '/Users/paulabetriu/Desktop/GitHub/SPICElib'; % library of
% SPICE functions (specially relevant for SPICE initialization in this
% program)
metricspath = '/Users/paulabetriu/Desktop/GitHub/science-opportunity/queries/geometric';

% Add paths
addpath(kernelpath);
addpath(resslib);
addpath(metricspath);

% Enable the access to the folders with the mosaic algorithms
addpath(genpath('/Users/paulabetriu/Desktop/GitHub/area-coverage-planning/'));

%% Case study
% Define mission and spacecraft (SPICE ID of the spacecraft)
mission= 'GALILEO';
sc     = 'GALILEO ORBITER';

% Choose instrument (SPICE ID of the instrument)
inst   = 'GLL_SSI';

% Define the planetary body (SPICE ID)
target = 'EUROPA';

% Planetary body modelization (DSK or ELLIPSOID)
method = 'ELLIPSOID';

% SPICE initialization with the relevant mission kernels
initMICE();
run(fullfile(kernelpath, lower(mission), 'inputkernels.m'));
loadKernels(fullK(METAKR));
%cspice_ldpool('/Users/paulabetriu/Desktop/MASTER/MUEA/TFM/SPICE/MICE/mice/kernels/gll36001.ti')

%% Definition of ROIs
% Pre-allocation of variables... 
stoptime = cspice_str2et('1998 MAY 30 00:00:00.000 TDB'); % mosaic end (max)
tobs  = 8.5; % [s] between observations
olapx = 20; % [%] of overlap in x direction
olapy = 20; % [%] of overlap in y direction
videosave = 0; % =1 saves a video with the evolution of the coverage map
% along the observation plan
count = 0;

% Regions of interest
% Pwyll Crater
count = count + 1;
roi = [ 78   -18;
        90   -14;
       109   -13;
       101   -23;
        99   -30;
        85   -33;
        71   -28;];
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).inittime = cspice_str2et('1998 MAR 29 12:38:00.000 TDB'); % closest approach
roistruct(count).name = "Pwyll Crater";

% Annwn Regio [lon, lat] = [40, 20]º
count = count + 1;
roi = [60 30;
       60 10;
       40 10;
       40 30;];
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).name = "Annwn Regio";
roistruct(count).inittime = cspice_str2et('1998 MAR 29 12:49:00.000 TDB'); % closest approach

% Niamh
count = count + 1;
roi = [150  25;
       150  15;
       135  15;
       135  25;]; % roi of roi polygon
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).name = "Niamh";
roistruct(count).inittime = cspice_str2et('1998 MAR 29 13:29:00.000 TDB');

% Cilix crater [lon, lat] = [180, 0]º;
count = count + 1;
roi = [-177  3;
       -177 -3;
        177 -3;
        177  3;];
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).name = "Cilix Crater";
roistruct(count).inittime = cspice_str2et('1998 MAR 29 13:40:00.000 TDB'); % closest approach

% Tara Regio
count = count + 1;
roi = [-55   20;
       -95   20;
       -95  -20;
       -55  -20;];
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).name = "Tara Regio";
roistruct(count).inittime = cspice_str2et('1998 MAR 29 14:00:00.000 TDB');

% Taliesin
count = count + 1;
roi = [-160  -10;
       -160  -30;
       -130  -30;
       -130  -10;]; % roi of roi polygon
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).name = "Taliesin";
roistruct(count).inittime = cspice_str2et('1998 MAR 29 14:21:00.000 TDB');
