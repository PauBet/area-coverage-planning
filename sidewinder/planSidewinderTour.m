function tour = planSidewinderTour(target, sc, t, roi, fprint0, ...
    gamma, olapx, olapy)
% This function discretizes and computes the coverage path of a certain
% ROI in order to build a mosaic image, adapted from [1]. Note that the
% coverage path is oriented as the width direction of the footprint
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        tour = planSidewinderTour(closestSide, roi, fprint0, gamma)
%
% Inputs:
%   > target:       SPICE string name of the target body
%   > sc:           SPICE string name of the spacecraft
%   > t:            time in TDB seconds past J2000 epoch
%   > roi:          matrix containing the roi of the ROI polygon. The
%                   vertex points are expressed in 2D. 
%       # roi(:,1) correspond to the x values of the roi
%       # roi(:,2) correspond to the y values of the roi
%   > fprint0:      struct containing the footprint parameters that are
%                   going to be used to define the grid discretization
%   > gamma:        origin of the grid, in latitudinal coordinates, in deg
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in deg
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in deg
% 
% Outputs:
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables
tour = {}; % list of planned observations
[~, targetFrame, ~] = cspice_cnmfrm(target); % target body-fixed reference 
% frame

% Calculate how is the spacecraft ground track position moving along the
% map. The coverage path does not only depend on the spacecraft position
% itself but also its velocity direction
subobs = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t,...
    targetFrame, 'NONE', sc);
[~, sclon, sclat] = cspice_reclat(subobs); % latitudinal coordinates
sclon = sclon*cspice_dpr; sclat = sclat*cspice_dpr; % [rad] to [deg]
subobs_ = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t + 60,...
    targetFrame, 'NONE', sc);
[~, sclon_, sclat_] = cspice_reclat(subobs_); % latitudinal coordinates
sclon_ = sclon_*cspice_dpr; sclat_ = sclat_*cspice_dpr; % [rad] to [deg]

% 2D grid discretization
grid = grid2D(fprint0, olapx, olapy, gamma, roi);

% Closest polygon side to the spacecraft's ground track position
cside = closestSide(target, sc, t, roi);

% The origin of the coverage path depends on the spacecraft ground track
% position
switch cside
    case {'up','down'} % Horizontal sweep
        if (sclon - sclon_) >= 0 
            bearing = false; % if the spacecraft is moving left (in the
            % topography map) then the coverage path should start at the
            % position furthest to the right
        else % right
            bearing = true; % if the spacecraft is moving right (in the
            % topography map) then the coverage path should start at the
            % position furthest to the left
        end
        for i=1:size(grid,1)
            % Sweep across latitude
            if isequal(cside, 'up')
                irow = i;
            else
                irow = size(grid, 1) - i + 1;
            end
            for j=1:size(grid, 2)
                if ~bearing
                    icol = size(grid, 2) - j + 1;
                else
                    icol = j;
                end
                if ~isempty(grid{irow, icol})
                    y = grid{irow, icol}(2);
                    x = grid{irow, icol}(1);
                    tour{end + 1} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each 
            % row sweeping, i.e. left (highest lon) to right (lowest lon) 
            % or vice versa
        end
    case {'right','left'} % Vertical sweep
        if (sclat - sclat_) < 0 % down
            bearing = false; % if the spacecraft is moving down (in the
            % topography map) then the coverage path should start at the
            % position furthest to the top
        else % up
            bearing = true; % if the spacecraft is moving up (in the
            % topography map) then the coverage path should start at the
            % position furthest to the bottom
        end
        for i=1:size(grid,2)
            % Sweep across longitude
            if isequal(cside, 'right')
                icol = size(grid, 2) - i + 1;
            else
                icol = i;
            end
            for j=1:size(grid, 2)
                if ~bearing
                    irow = j;
                else
                    irow = size(grid, 1) + 1 - j;
                end
                if ~isempty(grid{irow, icol})
                    y = grid{irow, icol}(2);
                    x = grid{irow, icol}(1);
                    tour{end + 1} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each
            % column sweeping, i.e. up (highest lat) to down (lowest
            % lat) or vice versa
        end
end

end