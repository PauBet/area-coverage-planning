% Mosaic comparison between the different heuristics

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
addpath(genpath(pwd));
addpath(genpath('/Users/paulabetriu/Desktop/GitHub/area-coverage-planning/'))

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
tcadence = 15; % [s] between observations
olapx = 20; % [%] of overlap in x direction
olapy = 20; % [%] of overlap in y direction
speedUp = 0;
count = 0;

% Regions of interest
count = count + 1;
roi = [25  25;
       25  -5;
       -5  -5;
       -5  25];
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).name = "ROI";
roistruct(count).inittime = cspice_str2et('1998 MAR 29 12:35:00.000 TDB'); % closest approach