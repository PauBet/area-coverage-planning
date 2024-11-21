clc; clear all;
% Revision of Sidewinder:
% Grid is going to be built in the camera frame, instead of the body-fixed
% frame. This is somewhat more complicated (it requires more calculations)
% but it could correct the spatial aberration that we presently see
% Re-implementation of sidewinder function with these new feature (and
% other code improvements)

% Load mission info (kernels, SPICE ids, etc.)
input_data_fig4;
roiname = char(lower(roistruct(1).name));
roiname(isspace(roiname)) = [];
name = ['post_process_',roiname];

% Define program iteration info
td = cspice_str2et('1998 MAR 29 12:00:00.000 TDB'); % initial observation ANNWN REGIO
%td = cspice_str2et('1998 MAR 29 13:40:00.000 TDB'); % initial observation TARA REGIO

% time
N = 150; % number of iterations
mkspan = zeros(1, N); % initialize mkspan array
coverage = zeros(1, N); % initialize coverage array
overlap  = zeros(1, N); % initialize coverage array
nfp = zeros(1, N); % initialize number of acquisitions array
step = 30; % time step in [sec]
et = linspace(1, N);

for i=1:N
    et(i) = td;

    % Grid Nibbler
    [A, fplist] = neighbour_placement_2(td, tcadence, inst, sc, ...
    target, roi, olapx, olapy, 3*1e-3);

    % % Plot tour
    % plotTour(A, fplist, roistruct, sc, target)
    % title(roistruct(1).name)
    % leg = legend('NumColumns', 2, 'Location', 'north');
    % leg.String(end) = [];
    
    % Get coverage, number of acquisitions and makespan
    [coverage(i), overlap(i)] = roicoverage(target, roi, fplist);
    mkspan(i) = fplist(end).t + tcadence - fplist(1).t;
    nfp(i) = length(fplist);
    
    % Next iteration
    td = td + step;    
end