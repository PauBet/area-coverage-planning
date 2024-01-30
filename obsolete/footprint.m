function fp = footprint(varargin)
% Given the spacecraft trajectory, this function computes the FOV
% projection onto the body surface, i.e. the footprint. The spacecraft
% orientation is determined by the lon, lat and theta angles.
% Assumption: the FOV is rectangular
% Note: this function is contained in the framework of the automated
% scheduler that optimizes the observation plan of a mission according to
% a set of scientific objectives. The spacecraft is assumed to be 3-axis 
% steerable and, therefore, CK kernels of the spacecraft (past missions)
% are avoided.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
%
% Usage:        fp = footprint(lon, lat, t, inst, sc, target, theta)
%
% Inputs:
%   > lon:      longitude coordinate of the target body at which the 
%               instrument boresight is pointing, in [deg]
%   > lat:      latitude coordinate of the target body at which the 
%               instrument boresight is pointing, in [deg]
%   > t:        time epoch in TDB seconds past J2000 epoch
%   > inst:     string SPICE name of the instrument
%   > sc:       string SPICE name of the spacecraft
%   > target:   string SPICE name of the target body
%   > theta:    roll angle of the instrument over its boresight axis, in
%               [deg]
%
% Outputs:
%   > fp:       struct containing main parameters of the footprint. In this
%               case, only the field 'vertices' is necessary.
%       # inst:          string SPICE name of the instrument
%       # sc:            string SPICE name of the spacecraft
%       # target:        string SPICE name of the target body
%       # t:             time epoch in TDB seconds past J2000 epoch
%       # bvertices:     matrix (N, 2) containing the N boundary vertices
%                        of the footprint polygon, in latitudinal 
%                        coordinates, in [deg].
%               + fp may contain more than one polygon, split because the
%                 actual footprint intersects the anti-meridian of the body
%                 surface. In this case, the polygons are separated in
%                 'vertices' by [NaN, NaN] Example:
%                         area1 = [350 30; 360 30; 360  0; 350 0];
%                         area2 = [ 0 30; 10 30; 10  0; 0  0];
%                         area = [area1; [NaN NaN]; area2];
%                         figure plot(polyshape(area(:,1), area(:,2));
%       # olon:          longitude coordinate of the target body at which
%                        the instrument boresight is pointing, in [deg]
%       # olat:          latitude coordinate of the target body at which 
%                        the instrument boresight is pointing, in [deg]
%       # sizex:         footprint size in one direction, in [deg]
%       # sizey:         footprint size in the other direction, in [deg]
%       # fovbsight:     FOV boresight in the body-fixed reference frame
%                        centered at the spacecraft position at time t
%       # fovbounds:     FOV bounds in the body-fixed reference frame
%                        centered at the spacecraft position at time t
%       # poleIntercept: string that determines if the footprint intercepts
%                        the north (='north pole') or south (='south pole')
%                        poles.
% 

%% Pre-allocate variables
lon     = varargin{1};
lat     = varargin{2};
t       = varargin{3};
inst    = varargin{4};
sc      = varargin{5};
target  = varargin{6};
theta   = varargin{7};
methodres  = 'highres';
pllflag = 0; % parallel execution
kernelPool = [];
if nargin > 7
    methodres = varargin{8};
%    pllflag = varargin{9};
%    kernelPool = varargin{10};
end

% If parallel execution is requested, the kernel pool needs to be
% initialized for each worker
if pllflag, initSPICE; end

% Footprint relevant fields
fp.inst      = inst;
fp.sc        = sc;
fp.target    = target;
fp.t         = t;
fp.bvertices = []; % footprint boundary vertices, in latitudinal 
% coordinates, in deg
fp.olon      = lon; % longitude value of FOV boresight projection onto the 
% body surface, in [deg]
fp.olat      = lat; % latitude value of FOV boresight projection onto the 
% body surface, in [deg]
fp.sizex     = 0; % horizontal size (longitude) [deg]
fp.sizey     = 0; % vertical size (latitude) [deg]
fp.fovbsight = []; % FOV boresight in target frame centered at the 
% spacecraft position
fp.fovbounds = []; % FOV bounds in target frame centered at the spacecraft
% position
fp.poleIntercept = []; % determines if the footprint intercepts the north
% pole (='north pole') or the south pole (='south pole')
fp.limb      = false; % boolean that defines if the FOV projects onto the
% planetary body's limb

%% Parameters
method = 'ELLIPSOID'; % assumption: ray intercept function is going to 
% model the target body as a tri-axial ellipsoid
lon = lon*cspice_rpd; % longitude from [deg] to [rad]
lat = lat*cspice_rpd; % latitude from [deg] to [rad]
theta = theta*cspice_rpd; % yaw angle from [deg] to [rad]
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % instrument FOV's boundary
    % vectors in the instrument frame
minx = min(bounds(1,:)); % minimum x focal plane
maxx = max(bounds(1,:)); % maximum x focal plane
miny = min(bounds(2,:)); % minimum y focal plane
maxy = max(bounds(2,:)); % maximum y focal plane
z = bounds(3,1); % z-coordinate of the boundary vectors
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
abcorr = 'LT'; % one-way light time aberration correction parameter.
% See cspice_spkpos for further information.
recpoint = cspice_srfrec(cspice_bodn2c(target), lon, lat); % rectangular
% coordinates of the target point in the body-fixed reference frame
instpos  = cspice_spkpos(sc, t, targetframe, abcorr, target); % rectangular
% coordinates of the instrument in the body-fixed reference frame

% Footprint vertices resolution
if isequal(methodres, 'lowres')
    N = 10;
else
    N = 200;
end

% Calculate instrument fov bounds (when this is pointing at the target 
% area)
%[bounds, boresight, pointingRotation, found] = instorient(inst, ...
%    target, lon*cspice_dpr, lat*cspice_dpr, sc, t, theta);
[bounds, boresight, pointingRotation, found] = instpointing(inst, ...
    target, lon*cspice_dpr, lat*cspice_dpr, sc, t, theta, 0);
if ~found
    return; % the point is not visible from the instrument's FOV, therefore
    % the function is exited and the footprint is returned empty
end

%% Field-of-view projection onto the target body
boundPoints = zeros(3,length(bounds)); % intercept points of the FOV
% boundary, in Cartesian coordinates, and in the body-fixed reference frame
% In its simplest form, the footprint should have the same number of
% vertices as boundaries has the instrument's FOV
irr = false; % boolean to define if the instrument FOV projection is not
% enclosed in the body surface, i.e. at least one of the boundary vectors
% do not intercept the body surface
intsec = false;
for i=1:length(bounds)
    [boundPoints(:,i), ~, ~, found] = cspice_sincpt(method, target, t,...
        targetframe, abcorr, sc, targetframe, bounds(:, i));
    % If the FOV boundary does not intercept the target...
    if ~found
        irr = true;
        fp.limb = true;
    else
        intsec = true;
    end
end

if ~intsec
    % If none of the 4 FOV boundaries intercept the target, and the target is
    % within FOV, then we're seeing the limb (might be easier to compute)
    % Parameters for cspice_limbpt function
    lbmethod = 'TANGENT/ELLIPSOID';
    [~, targetframe, ~] = cspice_cnmfrm(target); % body-fixed frame
    abcorr = 'XLT+S';
    corloc = 'CENTER';
    refvec = [0; 0; 1]; % first of the sequence of cutting half-planes
    ncuts  = 1e3; % number of cutting half-planes
    delrol = cspice_twopi() / ncuts; % angular step by which to roll the
    % cutting half-planes about the observer-target vector
    schstp = 1.0d-4; % search angular step size
    soltol = 1.0d-7; % solution convergence tolerance

    % Limb calculation with cspice_limbpt function
    [~, limb, ~, ~] = cspice_limbpt(lbmethod, target, t, targetframe, abcorr, ...
        corloc, sc, refvec, delrol, ncuts, schstp, soltol, ncuts); % limb
    % points expressed in targetframe ref frame
    surfPoints = limb';
else

    surfPoints = []; % matrix that saves the rectangular coordinates of the
    % intercept points between the FOV perimeter and the body surface
    count = 0; % counter of found intercept points
    if irr
        % Warning message
        if isequal(method, 'lowres')
            warning(['It is likely that the footprint contains the limb, ' ...
                'low resolution method may lead to significant inaccuracies'])
        end

        % Initialize variables
        maxfx = minx; maxfy = miny;
        minfx = maxx; minfy = maxy;

        % For those cases where the FOV does not completely contain the target,
        % a more refined search is going to be performed in order to define the
        % limits of the footprint
        for i=0:N
            % Vertical sweep of the focal plane
            x = (maxx - minx)*i/N + minx;
            for j=0:N
                y = (maxy - miny)*j/N + miny;
                vec = [x, y, z]';
                vec = pointingRotation*vec; % transform vector coordinates to
                % target frame
                corloc = 'SURFACE POINT'; % since alt is close to 0, there
                % should not be a significant difference between the target and
                % surface point correction locus (see cspice_tangpt)
                found = false; % found intercept
                [~, alt, ~, aux, ~, ~] = cspice_tangpt(method, target, t,...
                    targetframe, abcorr, corloc, sc, targetframe, vec);
                if alt < 10
                    % When the footprint contains the limb, its intercept is
                    % irregular, meaning that the boundary is not a smooth
                    % curve (the limb) but a set of scattered points with
                    % certain deviation around the limb. Besides this, since
                    % the research consists of a set of discretized points, we
                    % may not always find the limb intercept point. This
                    % depends on the resolution of the discretized refined
                    % mesh. To avoid incurring in excessive computational
                    % demands, instead of calculating the intercept point by
                    % refining the mesh, we calculate the tangent point. This
                    % is the closest surface point of the surface to the
                    % "intercepting" ray. When the ray actually intercepts the
                    % surface, the parameter 'alt', which is the distance
                    % between the tangent points and the surface, is equal to
                    % 0. We may find the limb "intercept" by finding those
                    % points where 'alt' is close or equal to 0.
                    % Future work: to be more precise.. try to minimize alt
                    % along the vertical sweep line
                    found = true;
                end

                if found && (y == miny || y == maxy || x == minx || x == maxx)
                    % if the vector intercepts the surface and is at the focal
                    % plane boundary...
                    count = count + 1;
                    surfPoints(count, :) = aux;

                    if y == miny
                        minfy = miny;
                    elseif y == maxy
                        maxfy = maxy;
                    end

                    if x == minx
                        minfx = minx;
                    elseif x == maxx
                        maxfx = maxx;
                    end
                elseif j>0 && found ~= old_found
                    % if the vector intercept status changes from the previous
                    % one, we're sweeping across the object's limb
                    count = count + 1;
                    if old_found
                        surfPoints(count, :) = old_surfPoint; % save the
                        % previous intercept
                    else
                        surfPoints(count, :) = aux; % save the current
                        % intercept
                    end

                    if x < minfx
                        minfx = x;
                    end
                    if y < minfy
                        minfy = y;
                    end
                    if x > maxfx
                        maxfx = x;
                    end
                    if y > maxfy
                        maxfy = y;
                    end
                end

                old_found = found; % save element intercept status
                old_surfPoint = aux; % save element intercept point
            end
        end
    else
        % The FOV projection is enclosed in the target surface
        % Close polygon
        boundPoints(:, end+1) = boundPoints(:, 1);
        % high resolution:
        for i=1:length(boundPoints)-1
            % linear (approximation) interpolation between vertices to define
            % the boundary of the footprint
            v = boundPoints(:, i+1) - boundPoints(:, i);
            lambda = linspace(0, 1, N); % line parametrization
            for l=1:N
                count = count + 1;
                surfPoints(count, 1) = boundPoints(1, i) + v(1)*lambda(l); % x
                surfPoints(count, 2) = boundPoints(2, i) + v(2)*lambda(l); % y
                surfPoints(count, 3) = boundPoints(3, i) + v(3)*lambda(l); % z
            end
        end
    end
end

if isempty(surfPoints), return; end % the FOV does not intercept with the
% object at any point of its focal plane

%% Conversion from rectangular to latitudinal coordinates of the polygon vertices
vertices = zeros(length(surfPoints), 2); % matrix that saves the
% latitudinal coordinates of the intercept points between the FOV
% perimeter and the body surface
[surfPoints(:,1), surfPoints(:,2), surfPoints(:,3)] = ...
    sortcw(surfPoints(:,1), surfPoints(:,2), surfPoints(:,3)); % sort
% polygon boundary vertices in clockwise order (for representation)
for i=1:length(surfPoints)
    [~, auxlon, auxlat] = cspice_reclat(surfPoints(i,:)'); % rectangular
    % to latitudinal coordinates
    vertices(i, 1) = auxlon*cspice_dpr(); % longitude in [deg]
    vertices(i, 2) = auxlat*cspice_dpr(); % latitude in [deg]
end
% Future work: surfPoints does not need to be saved, we could convert
% from rectangular to latitudinal inside the first loop, instead of
% doing separately. The reason why it is not is because we need to sort
% the vertices in clockwise order, and the sortcw algorithm for 2D does
% not work with non-convex polygons...

%% Check if the footprint contains north or south poles
npoleint = 0; % north pole intercept (=1 if the north pole is inside the
% instrument's FOV)
spoleint = 0; % south pole intercept (=1 if the south pole is inside the
% instrument's FOV)
bbounds = zeros(3, length(bounds));
for i=1:length(bounds)
    bbounds(:,i) = bounds(:,i) + instpos; % focal plane boundary points
end
p1 = bbounds(:,1);
p2 = bbounds(:,2);
p3 = bbounds(:,3);
pnormal = normalize(cross(p1 - p2, p1 - p3), 'norm');
plane = cspice_nvp2pl(pnormal, p1); % focal plane
plinst = [maxx maxy;
    maxx miny;
    minx miny;
    minx maxy];

% North pole check
rnpole = cspice_srfrec(cspice_bodn2c(target), 0, pi/2); % ray origin
v1 = instpos - rnpole; % ray direction
[spoint, ~, ~, found] = cspice_sincpt(method, target, t, targetframe,...
    abcorr, sc, targetframe, -v1); % check if the north pole is visible as
% seen from the instrument (the ray coming from the north pole could
% intersect the focal plane but actually be in the dark side of the body as
% seen from the spacecraft)
if found && norm(spoint - rnpole) < 1
    [int, pint] = cspice_inrypl(rnpole, v1, plane); % ray-plane intercept
    if int
        pint = pointingRotation\(pint - instpos); % intercept in the
        % instrument frame coordinates
        if inpolygon(pint(1), pint(2), plinst(:,1), plinst(:,2)) % find if
            % the intercept is contained in the focal plane boundaries
            npoleint = 1;
        end
    end
    if npoleint; fp.poleIntercept = 'north pole'; end
end

% South pole check
rspole = cspice_srfrec(cspice_bodn2c(target), 0, -pi/2);
v2 = instpos - rspole;
[spoint, ~, ~, found] = cspice_sincpt(method, target, t, targetframe,...
    abcorr, sc, targetframe, -v2); % check if the south pole is visible as
% seen from the instrument
if found && norm(spoint - rspole) < 1
    [int, pint] = cspice_inrypl(rspole, v2, plane); % ray-plane intercept
    if int
        pint = pointingRotation\(pint - instpos); % intercept in the
        % instrument frame coordinates
        if inpolygon(pint(1), pint(2), plinst(:,1), plinst(:,2)) % find if
            % the intercept is contained in the focal plane boundaries
            spoleint = 1;
        end
    end
    if spoleint; fp.poleIntercept = 'south pole'; end
end

%% Check if the footprint intersects the anti-meridian
% To ease the footprint representation on the topography map, we must
% consider the case where the footprint intercepts with the anti-meridian.
% If it does, we split the footprint in two polygons, cleaved by the line
% that the original footprint is crossing (a.m.)
[fp.bvertices(:, 1), fp.bvertices(:, 2)] = amsplit(vertices(:,1),...
    vertices(:,2)); % save footprint vertices

%% Save outputs
% Save fov parameters
fp.fovbounds = bounds; % save fov bounds in the body-fixed reference frame
fp.boresight = boresight; % save fov boresight vector in
% the body-fixed reference frame

% Calculate and save footprint size
xaxis = [maxx, 0,    1]';
yaxis = [0,    maxy, 1]';
xaxis = pointingRotation*xaxis;
yaxis = pointingRotation*yaxis;

if ~intsec
    fp.angle = 0;
    fp.sizex = 180;
    fp.sizey = 180;
    [clon, clat] = centroid(polyshape(fp.bvertices(:, 1), ...
        fp.bvertices(:, 2)));
    fp.clon = clon; fp.clat = clat;
    return;
else
    if ~irr
        [~, ~, ~, xpoint, ~, ~] = cspice_tangpt(method, target, t,...
            targetframe, abcorr, 'SURFACE POINT', sc, targetframe, xaxis);
        [~, ~, ~, ypoint, ~, ~] = cspice_tangpt(method, target, t,...
            targetframe, abcorr, 'SURFACE POINT', sc, targetframe, yaxis);

        % If the footprint does not contain the limb, then the centroid of the
        % footprint and the boresight projection are the same
        fp.clon = fp.olon; fp.clat = fp.olat;
    else
        % If the footprint contains the limb, the projection of the x/y axis,
        % if these do not intersect the body, will fall in the limb instead of
        % the midpoints of the footprint. Furthermore, the gamma (surface point
        % where the camera boresight is pointing at) will not necessarily
        % coincide with the centroid of the footprint.

        % We "adjust" the focal plane to the box that encloses the part of the
        % focal plane that actually contains the body and re-define the x/y
        % axes within that box
        deltafy = maxfy - minfy;
        deltafx = maxfx - minfx;
        xaxis = [maxfx,               miny + deltafy/2,    1]';
        yaxis = [minx + deltafx/2 ,   maxfy,               1]';
        xaxis = pointingRotation*xaxis;
        yaxis = pointingRotation*yaxis;

        [~, ~, ~, xpoint, ~, ~] = cspice_tangpt(method, target, t,...
            targetframe, abcorr, corloc, sc, targetframe, xaxis);
        [~, ~, ~, ypoint, ~, ~] = cspice_tangpt(method, target, t,...
            targetframe, abcorr, corloc, sc, targetframe, yaxis);

        % To compute the footprint's angle, we need to re-define the point
        % "recpoint", which should coincide with the footprint's centroid
        [lon, lat] = centroid(polyshape(fp.bvertices(:, 1), ...
            fp.bvertices(:, 2))); % footprint centroid -> new camera's
        % boresight pointing
        % Save centroid values
        fp.clon = lon; fp.clat = lat;

        lon = lon*cspice_rpd; lat = lat*cspice_rpd;
        recpoint = cspice_srfrec(cspice_bodn2c(target), lon, ...
            lat); % rectangular coordinates of the centroid

    end
end

% Compute the footprint angle between the longitude axis (projection) 
% and the footprint's x-axis
[~, lonx, latx] = cspice_reclat(xpoint);
xedge = [lonx - lon, latx - lat];
xedge = xedge./norm(xedge);
fp.angle = atan2d(xedge(2), xedge(1));

% Compute the footprint angle between the longitude axis (projection) 
% and the footprint's y-axis
[~, lony, laty] = cspice_reclat(ypoint);
yedge = [lony - lon, laty - lat];
yedge = yedge./norm(yedge);
angley = atan2d(yedge(2), yedge(1));
theta = abs(angley - fp.angle) - 90;

% Compute the footprint size as the angle separation between the
% boresight and the focal plane's x-axis projection in rectangular
% coordinates
if ~irr
    fp.sizex = 2*cspice_vsep(recpoint, xpoint)*cspice_dpr;
    fp.sizey = 2*cspice_vsep(recpoint, ypoint)*cspice_dpr;
else
    % Compute footprint size from the smallestBoundingBox (the method above
    % does not work well with limb projections)
    bbox = smallestBoundingBox(fp.bvertices(:, 1), fp.bvertices(:, 2), ...
        fp.angle);
    fp.sizex = bbox.size1;
    fp.sizey = bbox.size2;
end

fp.sizey = fp.sizey*abs(cosd(theta));

end