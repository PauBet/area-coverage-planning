function initMICE()

% get path to mice
HOMESPICE = getHomeSpice;

% Add mice and lib folders from mice
addpath(fullfile(HOMESPICE,'src/mice/'))
addpath(fullfile(HOMESPICE,'lib'))

end