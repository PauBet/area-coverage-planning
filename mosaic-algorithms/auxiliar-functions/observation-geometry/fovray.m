function visible = fovray(inst, target, obs, et, lon, lat, varargin)
% Pre-allocate variables
[~, rframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE

if nargin > 6
    % Determine spacecraft attitude according to lon, lat boresight
    % pointing
    blon = varargin{1};
    blat = varargin{2};

    % Project latitudinal coordinates on body surface to focal plane
    instpoint = topo2inst([lon, lat], blon, blat, target, obs, inst, et);
    if isempty(instpoint)
        % Point is not visible
        visible = false;
        return;
    end
    
    % If point is visible, let's check if it is inside the FOV limits
    % Retrieve FOV parameters
    [~, ~, ~, bounds] = ...
        cspice_getfov(cspice_bodn2c(inst), 4); % instrument FOV's boundary
    % vectors in the instrument frame
    % Get min-max FOV boundaries in the focal plane
    maxx = max(bounds(1, :)); minx = min(bounds(1, :));
    maxy = max(bounds(2, :)); miny = min(bounds(2, :));

    if instpoint(1) >= minx && instpoint(1) <= maxx && instpoint(2) >= miny ...
            && instpoint(2) <= maxy
        % Point is visible and in FOV
        visible = true;
    else
        % Point is visible but not in FOV
        visible = false;
    end
else
    % Check point in FOV by retrieving spacecraft attitude from CK
    abcorr = 'NONE';
    raydir = -trgobsvec([lon lat], et, target, sc);
    visible = cspice_fovray(inst, raydir, rframe, abcorr, obs, et);
end

end