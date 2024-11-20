clc; close all; clear all;
% Calculation of coverage and makespan for the same ROI over time with
% Online Frontier

% Load mission info (kernels, SPICE ids, etc.)
input_data_fig4;
mosaic = 'onlinefrontier';
roiname = char(lower(roistruct(1).name));
roiname(isspace(roiname)) = [];
name = ['post_process_',roiname];

% Define program iteration info
td = cspice_str2et('1998 MAR 29 12:10:00.000 TDB'); % initial observation
% time
N = 15; % number of iterations
mkspan = zeros(1, N); % initialize mkspan array
coverage = zeros(1, N); % initialize coverage array
nfp = zeros(1, N); % initialize number of acquisitions array
step = 30; % time step in [sec]

for i=1:N
    % Online Frontier
    [A, fplist] = frontierRepair(td, stoptime, tcadence, inst, sc, ...
        target, roi, olapx, olapy, 3*1e-3);

    % % Plot tour
    % plotTour(A, fpList, roistruct, sc, target)
    % title(roistruct(1).name)
    % leg = legend('NumColumns', 2, 'Location', 'north');
    % leg.String(end) = [];
    
    % Get coverage, number of acquisitions and makespan
    coverage(i) = roicoverage(target, roi, fplist);
    mkspan(i) = fplist(end).t + tcadence - fplist(1).t;
    nfp(i) = length(fplist);
    
    % Next iteration
    td = td + step;    
end

figure
plot(1:N, mkspan, 'b+:')

figure
plot(1:N, coverage, 'r^:')

figure
plot(1:N, nfp, 'g*:')

% FOM post-process
%run(name);