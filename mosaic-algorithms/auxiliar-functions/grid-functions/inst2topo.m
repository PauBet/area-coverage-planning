function grid_topo = inst2topo(grid, lon, lat, target, sc, inst, et)
% This function transforms a set of points from the instrument frame to the
% topographic coordinate system (latitude and longitude on the target body)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        grid_topo = inst2topo(grid, lon, lat, target, sc, inst, et)
%
% Inputs:
%   > grid:         cell array of grid points in instrument frame coordinates.
%                   Each cell contains a 2D point [x, y] representing a 
%                   location in the instrument frame, or is empty if no 
%                   observation point is defined
%   > lon:          longitude of the observation point or area center, in
%                   [deg]
%   > lat:          latitude of the observation point or area center, in
%                   [deg]
%   > target:       string name of the target body
%   > sc:           string name of the spacecraft
%   > inst:         string name of the instrument
%   > et:           ephemeris time, TDB seconds past J2000 epoch
% 
% Outputs:
%   > grid_topo:    cell array of the input grid points transformed to 
%                   topographic coordinates on the target body. Each cell 
%                   contains a 2D point [lon, lat], in [deg]

% Pre-allocate
[~, ~, rotmat] = instpointing(inst, target, sc, et, lon, lat);
grid_topo = cell(size(grid));
method = 'ELLIPSOID';
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE

% Convert grid into topographical coordinates
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        instp = grid{i, j}; % retrieve current point in instrument frame
        if ~isempty(instp) && ~any(isnan(instp), 'all')
            p = zeros(3, 1); % initialize 3D points for conversion
            p(1:2) = instp; p(3) = 1; % assign [x, y] coordinates 
            p_body = rotmat*p; % set z to 1 for surface projection 
            % calculation

            % Compute surface intersection of the point on the target body
            [xpoint, ~, ~, found] = cspice_sincpt(method, target, et,...
                targetframe, 'NONE', sc, targetframe, p_body);
            if found
                % Convert rectangular coordinates to latitudinal
                [~, lon, lat] = cspice_reclat(xpoint);
                grid_topo{i, j} = [lon*cspice_dpr, lat*cspice_dpr];
            else
                disp("Point not visible from the instrument")
            end
        end
    end
end
end