function [vroi, poly1] = visibleroi(roi, et, target, obs)
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
method = 'TANGENT/ELLIPSOID';
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
[~, limb, ~, ~] = cspice_limbpt(method, target, et, targetframe, abcorr, ...
    corloc, obs, refvec, delrol, ncuts, schstp, soltol, ncuts); % limb
% points expressed in targetframe ref frame
[~, lblon, lblat] = cspice_reclat(limb); % conversion from
% rectangular to latitudinal coordinates
lblon = lblon*cspice_dpr;
lblat = lblat*cspice_dpr;

% Check a.m. split
ind2 = find(diff(sort(lblon)) >= 180, 1); % find the discontinuity 
if ~isempty(ind2)
    [lblon, lblat] = amsplit(lblon', lblat');
end

% roi and limb intersection
poly1 = polyshape(lblon, lblat);
poly2 = polyshape(roi(:, 1), roi(:, 2));
inter = intersect(poly1, poly2);

% output visible roi
vroi  = inter.Vertices;

end