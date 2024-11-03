function [vroi, inter, flag] = visibleroi(roi, et, target, obs)
% Given a target area (region-of-interest) on a planetary body surface and 
% an observer, this function calculates the portion of the former that is
% visible
% Note: if roi intercepts the anti-meridian line, this function returns 
% the visible roi polygon divided by this line (see amsplit function)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
%               Diego Andía  (UPC/ESEIAAT)
% Date:         06/2023
%
% Usage:        vroi = visibleroi(roi, et, target, obs)
%
% Inputs:
%   > roi:     matrix containing the vertices of the ROI polygon. The
%              vertex points are expressed in 2D, in latitudinal
%              coordinates [º]
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices
%   > et:      limb projection time, in seconds past J2000 epoch
%   > target:  string name of the target body (SPICE ID)
%   > obs:     string name of the observer (SPICE ID)
%
% Output:
%   > vroi:    matrix containing the vertices of the intersection between
%              the input ROI polygon and the limb projection, i.e., the 
%              visible portion of the ROI on the body surface as seen from 
%              the observer

% Previous anti-meridian intersection check...
ind1 = find(diff(sort(roi(:, 1))) >= 180, 1); % find the discontinuity index
if ~isempty(ind1)
    [aux(:, 1), aux(:, 2)] = amsplit(roi(:, 1), roi(:, 2));
    roi = aux;
end

% Parameters for cspice_limbpt function
flag = false; % assume the target area is visible from the instrument
method = 'TANGENT/ELLIPSOID';
[~, targetframe, ~] = cspice_cnmfrm(target); % body-fixed frame
abcorr = 'LT+S';
corloc = 'CENTER';
refvec = [0; 0; 1]; % first of the sequence of cutting half-planes
ncuts  = 1e3; % number of cutting half-planes
delrol = cspice_twopi() / ncuts; % angular step by which to roll the
% cutting half-planes about the observer-target vector
schstp = 1.0d-4; % search angular step size
soltol = 1.0d-7; % solution convergence tolerance

% Limb calculation with cspice_limbpt function
[~, limb, ~, ~] = cspice_limbpt(method, target, et, targetframe, abcorr, ...
    corloc, obs, refvec, delrol, ncuts, schstp, soltol, ncuts); % limb
% points expressed in targetframe ref frame
[~, lblon, lblat] = cspice_reclat(limb); % conversion from
% rectangular to latitudinal coordinates
lblon = lblon*cspice_dpr;
lblat = lblat*cspice_dpr;

% Check for north/south pole
northpole = false;
southpole = false;
% Check north-pole:
srfpoint = [0 90];
angle = emissionang(srfpoint, et, target, obs);
if angle < 90, northpole = true; end
% Check south-pole:
srfpoint = [0 -90];
angle = emissionang(srfpoint, et, target, obs);
if angle < 90, southpole = true; end

% Case 1.
if ~northpole && ~southpole
    % Check a.m. split
    ind2 = find(diff(sort(lblon)) >= 180, 1); % find the discontinuity
    if ~isempty(ind2)
        [lblon, lblat] = amsplit(lblon', lblat');
    end
    % Check if we are keeping the correct polygon (full disk polygons may be
    % misleading, we can only guarantee through emission angle check)
    exit = 0;
    while ~exit
        randPoint = [randi([-180 180]), randi([-90 90])];
        if inpolygon(randPoint(1), randPoint(2), lblon, lblat)
            angle = emissionang(randPoint, et, target, obs);
            if angle < 85
                exit = 1;
                poly1 = polyshape(lblon, lblat);
            end
        else
            angle = emissionang(randPoint, et, target, obs);
            if angle < 85
                exit = 1;
                % This calculation is approximated, we should find a better way
                % to find the complementary
                %% [Future work]
                lonmap = [-180 -180 180 180];
                latmap = [-90    90  90 -90];
                polymap = polyshape(lonmap, latmap);
                poly1 = polyshape(lblon, lblat);
                poly1 = subtract(polymap, poly1);
            end
        end
    end
else
    % Case 2.
    [lblon, indsort] = sort(lblon);
    lblat = lblat(indsort);
    if northpole || southpole
        % Include northpole to close polygon
        auxlon = lblon; auxlat = lblat;
        lblon  = zeros(1, length(auxlon) + 2); 
        lblat = zeros(1, length(auxlat) + 2);
        if northpole
            lblon(1) = -180; lblat(1) = 90;
            lblon(end) = 180; lblat(end) = 90;
        else
            lblon(1) = -180; lblat(1) = -90;
            lblon(end) = 180; lblat(end) = -90;
        end
        lblon(2:length(lblon)-1) = auxlon; lblat(2:length(lblat)-1) = auxlat;
    end
    poly1 = polyshape(lblon, lblat);
end

% roi and limb intersection
poly2 = polyshape(roi(:, 1), roi(:, 2));
inter = intersect(poly1, poly2);

% output visible roi
vroi  = inter.Vertices;

% visibility flag
if isempty(vroi), flag = true; end

end