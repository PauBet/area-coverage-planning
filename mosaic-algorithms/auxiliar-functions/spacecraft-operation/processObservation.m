function [A, tour, fpList, poly1, t, empty] = processObservation(A, tour, ...
    fpList, poly1, t, slewRate, tobs, amIntercept, inst, sc, target, ...
    resolution)
% This function handles the processing of an observation point by computing
% its footprint, updating the list of completed observations and adjusting 
% the remaining area to be covered. It also accounts for anti-meridian 
% interception and updates the time for the next observation based on slew 
% rate and observation time.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2023
%
% Usage:        [A, tour, fpList, poly1, t] = processObservation(A, tour, ...
%                   fpList, poly1, t, slewRate, tobs, amIntercept, inst, sc, target, ...
%                   resolution)
%
% Inputs:
%   > A:            cell matrix of successive instrument observations,
%                   sorted in chronological order.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > tour:         array of remaining observation points in the tour, in
%                   latitudinal coordinates [º]
%   > fpList:       list of footprint structures detailing the observation
%                   metadata and coverage
%   > poly1:        current polygon shape of the uncovered area on the 
%                   target body
%   > t:            current time in ephemeris seconds past J2000 epoch
%   > slewRate:     rate at which the spacecraft (or instrument platform)
%                   can slew between observations, in [º/s]
%   > tobs:         observation time, i.e. the minimum time that the 
%                   instrument needs to perform an observation, in seconds
%   > amIntercept:  boolean flag indicating if the anti-meridian is
%                   intercepted by the observation path
%   > inst:         string name of the instrument
%   > sc:           string name of the spacecraft
%   > target:       string name of the target body
%   > resolution:   string defining the resolution setting, affecting the
%                   footprint calculation. It could be either 'lowres' or
%                   'highres'. See footprint function for further
%                   information
%
% Outputs:
%   > A, tour, fpList, poly1, t: updated variables

% Previous check...
if isempty(tour), empty = true; return; end

% Compute the footprint of each point in the tour successively and
% subtract the corresponding area from the target polygon
a = tour{1}; % observation
tour(1) = []; % delete this observation from the planned tour
empty = false;

% Check a.m. intercept...
if a(1) > 180, a(1) = a(1) - 360; end

% Compute the observation's footprint
fprintf('Computing %s FOV projection on %s at %s...', inst, ...
    target, cspice_et2utc(t, 'C', 0));
fprinti = footprint(t, inst, sc, target, resolution, a(1), a(2), 0);
% Body-fixed to inertial frame

if ~isempty(fprinti.bvertices)
    fprintf('\n')
    % Check a.m. intercept
    if amIntercept
        aux = fprinti; ind = aux.bvertices(:, 1) < 0;
        aux.bvertices(ind, 1) = aux.bvertices(ind, 1) + 360;
        poly2 = polyshape(aux.bvertices);
    else
        poly2 = polyshape(fprinti.bvertices); % create footprint
        % polygon
    end
    
    % Check footprint-ROI intersect
    targetpshape = poly1;
    areaT = area(targetpshape);
    inter = subtract(targetpshape, poly2);
    areaI = area(inter);
    areaInter = areaT - areaI;
    fpArea = area(poly2);

    if areaInter/fpArea == 0
        empty = true;
    else
        A{end + 1} = a; % add it in the list of planned observations
        poly1 = subtract(poly1, poly2); % update uncovered area

        % Save footprint struct
        fpList(end + 1) = fprinti;

        % New time iteration
        if ~isempty(tour)
            p1 = [fprinti.olon, fprinti.olat];
            p2 = [tour{1}(1), tour{1}(2)];
            t = t + tobs + slewDur(p1, p2, t, tobs, inst, target, sc, slewRate);
        end
    end
else
    empty = true;
end

if empty
    fprintf(' Surface not reachable\n')
end

end