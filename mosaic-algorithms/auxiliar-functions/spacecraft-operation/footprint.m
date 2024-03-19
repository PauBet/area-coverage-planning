function fp = footprint(t, inst, sc, target, res, varargin)
% Given the spacecraft trajectory, this function computes the FOV
% projection onto the body surface, i.e. the footprint. The spacecraft
% orientation is determined by the lon, lat and theta angles.
% Assumption: the FOV is rectangular
% Note: this function is contained in the framework of the automated
% scheduler that optimizes the observation plan of a mission according to
% a set of scientific objectives. 
% [Warning]: limb projections on a topographical map might be very
% irregular. FUTURE WORK
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% Version:      2
% Last update:  12/2023
%
% Usage:        fp = footprint(t, inst, sc, target, res)
%               fp = footprint(t, inst, sc, target, res, lon, lat)
%
% Inputs:
%   > t:        time epoch in TDB seconds past J2000 epoch
%   > inst:     string SPICE name of the instrument
%   > sc:       string SPICE name of the spacecraft
%   > target:   string SPICE name of the target body
%   > lon:      longitude coordinate of the target body at which the
%               instrument boresight is pointing, in [deg]
%   > lat:      latitude coordinate of the target body at which the
%               instrument boresight is pointing, in [deg]
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
%       # fovbsight:     FOV boresight in the body-fixed reference frame
%                        centered at the spacecraft position at time t
%       # fovbounds:     FOV bounds in the body-fixed reference frame
%                        centered at the spacecraft position at time t
%

%% Pre-allocate variables

% Ray intercept function parameters
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
abcorr = 'LT+S'; % one-way light time aberration correction parameter
method = 'ELLIPSOID'; % assumption: ray intercept function is going to
% model the target body as a tri-axial ellipsoid
surfPoints = []; % matrix that saves the rectangular coordinates of the
% intercept points between the FOV perimeter and the body surface
count = 0; % counter of found intercept points
geom = false;

% Definition of instrument pointing (3-axis steerable or constrained)
if nargin == 5 % instrument pointing is provided by (and retrieve from) a ck
    ckpointing = true;
elseif nargin > 5 % instrument is considered 3-axis steerable
    lon = varargin{1};
    lat = varargin{2};
    ckpointing = false;
    if nargin == 8
        geom = varargin{3};
    end
end

% Definition of footprint resolution
if isequal(res, 'lowres') % footprint vertices resolution
    N = 10; % number of intercept search per side
elseif isequal(res, 'highres')
    N = 500;
else
    error("Invalid resolution method")
end

% Initialize footprint struct fields
fp.inst      = inst;
fp.sc        = sc;
fp.target    = target;
fp.t         = t;
fp.bvertices = []; % footprint boundary vertices, in latitudinal
% coordinates, in deg
fp.olon      = NaN; % longitude value of FOV boresight projection onto the
% body surface, in [deg]
fp.olat      = NaN; % latitude value of FOV boresight projection onto the
% body surface, in [deg]
fp.fovbsight = []; % FOV boresight in target frame centered at the
% spacecraft position
fp.fovbounds = []; % FOV bounds in target frame centered at the spacecraft
% position
fp.limb      = 'none'; % boolean that defines if the FOV projects onto the
% planetary body's limb
fp.recVertices = []; % matrix that contain the rectangular coordinates of
% the footprint boundary vertices
fp.angle     = NaN;
fp.width     = NaN;
fp.height    = NaN;

%% Calculate instrument orientation
if ckpointing % constrained pointing
    [bounds, boresight, pointingRotation, found, lon, lat] = instpointing(inst, ...
        target, sc, t);
else % 3-axis steerable
    [bounds, boresight, pointingRotation, found] = instpointing(inst, ...
        target, sc, t, lon, lat);
end
if ~found
    return; % the point is not visible from the instrument's FOV, therefore
    % the function is exited and the footprint is returned empty
end
fp.fovbounds = bounds; % save fov bounds in the body-fixed reference frame
fp.fovbsight = boresight; % save fov boresight vector in the body-fixed 
% reference frame

% Update boresight pointing (in latitudinal coordinates)
fp.olon = lon; % longitude value of FOV boresight projection onto the
% body surface, in [deg]
fp.olat = lat; % latitude value of FOV boresight projection onto the
% body surface, in [deg]

%% Project instrument FOV onto the body surface

% Retrieve FOV parameters
[~, ~, ~, bounds] = ...
    cspice_getfov(cspice_bodn2c(inst), 4); % instrument FOV's boundary
% vectors in the instrument frame
minx = min(bounds(1,:)); % minimum x focal plane
maxx = max(bounds(1,:)); % maximum x focal plane
miny = min(bounds(2,:)); % minimum y focal plane
maxy = max(bounds(2,:)); % maximum y focal plane
z = bounds(3,1); % z-coordinate of the boundary vectors

boundPoints = zeros(3,length(bounds)); % intercept points of the FOV
% boundary, in Cartesian coordinates, and in the body-fixed reference frame
% In its simplest form, the footprint should have the same number of
% vertices as boundaries has the instrument's FOV
fp.limb = 'none'; % boolean to define if the instrument FOV projection is not
% enclosed in the body surface, i.e. at least one of the boundary vectors
% do not intercept the body surface
intsec = false; % boolean that indicates if at least one of the boundary 
% vectors intercepts the body surface
for i=1:length(fp.fovbounds)
    [boundPoints(:,i), ~, ~, found] = cspice_sincpt(method, target, t,...
        targetframe, abcorr, sc, targetframe, fp.fovbounds(:, i));
    % If the FOV boundary does not intercept the target...
    if ~found
        % If at least one FOV boundary does not intercept the
        % object's surface, then we're seeing (at least, partially)
        % the limb
        fp.limb = 'partial';
    else
        intsec = true;
    end
end

% Assume that, if there are no intercepts of the FOV boundaries, the FOV 
% projection is likely to contain the total limb of the body
% In the next step, we will find out if this assumption is correct
if ~intsec
    fp.limb = 'total';
end

% When the footprint is likely to contain the limb, perform a more
% accurate search in order to conclude if the FOV intercepts the
% body at some point
if isequal(fp.limb, 'total'), refineFOVsearch(); end

%% Compute footprint
if isequal(fp.limb, 'none')
    % FOV projects entirely on the body surface
    inFOVprojection();
elseif isequal(fp.limb, 'partial')
    % FOV contains partially the limb
    plimbFOVprojection();
else
    % FOV contains the whole body (total limb)
    tlimbFOVprojection();
end

if isempty(surfPoints), return; end % the FOV does not intercept with the
% object at any point of its focal plane

% Save values
fp.recVertices = surfPoints;

%% Conversion from rectangular to latitudinal coordinates of the polygon vertices
% and geometry computation
footprint2map();

    function refineFOVsearch()
        % Perform a perimetral search of the FOV to find out if the FOV
        % contains totally or partially the body limb
        Nl = 20;
        found = false;
        for ii=0:1
            % Vertical sweep of the focal perimeter
            x = (maxx - minx)*ii + minx;
            for jj=0:Nl
                y = (maxy - miny)*jj/Nl + miny;
                vec = [x, y, z]';
                vec = pointingRotation*vec; % transform vector coordinates to
                % target frame
                [~, ~, ~, found] = cspice_sincpt(method, target, t,...
                    targetframe, abcorr, sc, targetframe, vec);
                if found
                    fp.limb = 'partial';
                    break;
                end
            end
            if found, break; end
        end

        % If intercept still has not been found...
        if ~found
            for ii=0:Nl
                % Horizontal sweep of the focal perimeter
                x = (maxx - minx)*ii/Nl + minx;
                for jj=0:1
                    y = (maxy - miny)*jj/Nl + miny;
                    vec = [x, y, z]';
                    vec = pointingRotation*vec; % transform vector coordinates to
                    % target frame
                    [~, ~, ~, found] = cspice_sincpt(method, target, t,...
                        targetframe, abcorr, sc, targetframe, vec);
                    if found
                        fp.limb = 'partial';
                        break;
                    end
                end
                if found, break; end
            end
        end
    end

    function inFOVprojection()
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

    function plimbFOVprojection()

        % Warning message
        if isequal(res, 'lowres')
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
                if alt < 15
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

    end

    function tlimbFOVprojection()
        % Compute limb with SPICE function (easier)
        % Parameters for cspice_limbpt function
        lbmethod = 'TANGENT/ELLIPSOID';
        abcorr = 'XLT+S';
        corloc = 'CENTER';
        refvec = [0; 0; 1]; % first of the sequence of cutting half-planes
        ncuts  = 2e3; % number of cutting half-planes
        delrol = cspice_twopi() / ncuts; % angular step by which to roll the
        % cutting half-planes about the observer-target vector
        schstp = 1.0d-6; % search angular step size
        soltol = 1.0d-10; % solution convergence tolerance

        % Limb calculation with cspice_limbpt function
        [~, limb, ~, ~] = cspice_limbpt(lbmethod, target, t, targetframe, abcorr, ...
            corloc, sc, refvec, delrol, ncuts, schstp, soltol, ncuts); % limb
        % points expressed in targetframe ref frame
        surfPoints = limb';
    end

    function footprint2map()

        % Pre-allocate variables
        vertices = zeros(length(surfPoints), 2); % matrix that saves the
        % latitudinal coordinates of the intercept points between the FOV
        % perimeter and the body surface

        % Sort points
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

        %% Check if the footprint intersects the anti-meridian
        % To ease the footprint representation on the topography map, we must
        % consider the case where the footprint intercepts with the anti-meridian.
        % If it does, we split the footprint in two polygons, cleaved by the line
        % that the original footprint is crossing (a.m.)
        [fp.bvertices(:, 1), fp.bvertices(:, 2)] = amsplit(vertices(:,1),...
            vertices(:,2)); % save footprint vertices
        
        if geom
            % Get minimum width direction and size
            [angle, width, height] = minimumWidthDirection(fp.bvertices(:, 1), ...
                fp.bvertices(:, 2));
            fp.angle = angle;
            fp.width  = width;
            fp.height = height;
        end

    end

end