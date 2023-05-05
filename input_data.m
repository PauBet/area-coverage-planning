%% ********** input data **********
% Paula Betriu - March 2023

%% Relevant paths
kernelpath = 'input';
resslib    = '/Users/paulabetriu/Desktop/GitHub/RESSlib'; % library of
% SPICE functions (specially relevant for SPICE initialization in this
% program)
metricspath = '/Users/paulabetriu/Desktop/GitHub/autonomous-scheduling/physical-quantities';

% Add paths
addpath(kernelpath);
addpath(resslib);
addpath(metricspath);

% Enable the access to the folders with the mosaic algorithms
addpath(genpath(pwd));

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
run(fullfile(kernelpath, lower(mission), 'inputkernels.m'));
initSPICEv(fullK(METAKR));
%cspice_ldpool('/Users/paulabetriu/Desktop/MASTER/MUEA/TFM/SPICE/MICE/mice/kernels/gll36001.ti')