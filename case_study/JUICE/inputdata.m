% Input data (from user) containing relevant information for the software

% Add necessary paths
addpath(genpath('/Users/paulabetriu/Desktop/GitHub/area-coverage-planning'));

% Pre-allocate variables
global figpath; %#ok
global pypath; %#ok
pypath = '/usr/local/bin/python3'; % python command
% Figures of merit are going to be saved in a folder named as the time of
% execution:
figpath = fullfile(pwd, 'figures-of-merit', ...
    string(datetime('now','Format','yyyyMMddHHmm')));
if ~isfolder(figpath), mkdir(fullfile(figpath)); end

%% Mission information
% Pre-allocate variables
obsinfo.target = 'GANYMEDE';
obsinfo.sc     = 'JUICE';
obsinfo.inst   = 'JUICE_JANUS';

% Initialize SPICE
initMICE()
% load metakernel
cspice_furnsh('/Users/paulabetriu/Desktop/phd/misc/juice-misc/juice/kernels/mk/juice_plan_local.tm')

%% Case study
% Define mission and spacecraft (SPICE ID of the spacecraft)
mission= 'JUICE';
sc     = 'JUICE';

% Choose instrument (SPICE ID of the instrument)
inst   = 'JUICE_JANUS';

% Define the planetary body (SPICE ID)
target = 'GANYMEDE';

% Planetary body modelization (DSK or ELLIPSOID)
method = 'ELLIPSOID';

%% Definition of ROIs
% Pre-allocation of variables... 
stoptime = cspice_str2et('2035 MAY 30 00:00:00.000 TDB'); % mosaic end (max)
tobs  = 8.5; % [s] between observations
olapx = 20; % [%] of overlap in x direction
olapy = 20; % [%] of overlap in y direction
videosave = 0; % =1 saves a video with the evolution of the coverage map
% along the observation plan
count = 0;

% Regions of interest
% Xibalba Sulcus
count = count + 1;
roi = [ -90   31;
        -68   31;
        -68   46;
        -90   46];
roistruct(count).vertices = roi;
[cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
roistruct(count).cpoint = [cx, cy];
roistruct(count).inittime = 1089532300.1838222; % closest approach
roistruct(count).name = "JUICE_ROI_GAN_1_0_02";