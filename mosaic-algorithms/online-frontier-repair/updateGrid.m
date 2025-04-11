function [seed, inst_grid, inst_tour, topo_tour] = updateGrid(roi, ...
    inst_tour, inst_grid, grid_dirx, grid_diry, cx, cy, olapx, olapy, ...
    insweepDir1, insweepDir2, seed, old_seed, gamma, et, inst, sc, target)
% This function dynamically updates the grid of observations by 
% incorporating new observation points, adjusting for changes in the 
% observation geometry, and removing points that no longer contribute 
% to covering the region of interest (ROI). It uses a boustrophedon pattern
% for traversal. Adapted from [1].
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2023
%
% Usage:        [seed, inst_grid, inst_tour, topo_tour] = updateGrid2(roi, ...
%                   inst_tour, inst_grid, grid_dirx, grid_diry, angle0, cx, cy, olapx, olapy, ...
%                   indir1, indir2, seed, old_seed, gamma, et, inst, sc, target)
%
% Inputs:
%   > roi:          matrix containing the vertices of the uncovered area
%                   of the ROI polygon. The vertex points are expressed in 
%                   2D, in latitudinal coordinates [º]
%   > inst_tour:    tour path in instrument frame coordinates
%   > inst_grid:    grid of potential observation points in instrument
%                   frame coordinates
%   > grid_dirx:    direction of the grid along the x-axis in the
%                   instrument frame
%   > grid_diry:    direction of the grid along the y-axis in the
%                   instrument frame
%   > cx, cy:       centroid coordinates of the roi
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in percentage (height)
%   > insweepDir1, insweepDir2: initial directions defining the sweep of the
%                           Boustrophedon decomposition, derived according 
%                           to the spacecraft's position with respect to
%                           the ROI
%   > seed:         current starting point for the grid update
%   > old_seed:     starting point from the previous iteration
%   > gamma:        current observation point
%   > et:           current time in ephemeris seconds past J2000 epoch
%   > inst:         string name of the instrument
%   > sc:           string name of the spacecraft
%   > target:       string name of the target body
%
% Outputs:
%   > seed, inst_grid, inst_tour: updated variables
%   > topo_tour:    tour path in topographical coordinates (lat/lon on the
%                   target body), in [deg]
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables...
persistent fpref;
persistent pointing0;
persistent sweepDir1;
persistent sweepDir2;
if isempty(sweepDir1)
    sweepDir1 = insweepDir1; sweepDir2 = insweepDir2;
end
if isempty(pointing0), [pointing0(1), pointing0(2)] = deal(cx, cy); end
cin = {}; % set of new observations (inside or outside from 'tour')
cind= {}; % map indices of the new observation (potential placement in the
% map)
cout = {}; % set of disposable observations (inside or outside from 'tour')
N    = {}; % set of new tiles (outside from 'tour')
Nind = {}; % map indices of the new tiles
X    = {}; % set of disposable tiles (inside of 'tour')
epsilon = 0.01;

% Build reference tile (it's always going to be the same in the subsequent
% calls)
if isempty(fpref)
    % We need to craft the tile reference that we are going to use 
    % throughout the heruristics operations. To avoid repetitions, we
    % declare it as a persistent variable that is only going to be
    % calculated in the first call
    [~, ~, ~, bounds] = ...
        cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
    maxx = max(bounds(1, :)); minx = min(bounds(1, :));
    maxy = max(bounds(2, :)); miny = min(bounds(2, :));
    width = maxx - minx; height = maxy - miny;
    fpref.width = width;
    fpref.height = height;
    xlimit = [-width/2 width/2]; ylimit = [-height/2 height/2];
    xbox = xlimit([1 1 2 2 1]); ybox = ylimit([1 2 2 1 1]);
    fpref.bvertices(:, 1) = xbox;
    fpref.bvertices(:, 2) = ybox;
end

% Project ROI topographical coordinates to instrument's focal plane
targetArea = topo2inst(roi, cx, cy, target, sc, inst, et);
targetpshape = polyshape(targetArea(:,1), targetArea(:,2)); % polyshape

% [Future work]: orientation angle may change over the course of the mosaic
%targetArea = topo2inst(roi, gamma_topo(1), gamma_topo(2), target, sc, ...
%    inst, et); % current roi coordinates in the instrument reference frame,
% when the instrument is pointing at the current grid origin point (next
% observation)
% Since we are not calculating the grid again, but rather updating a
% reference, we need to have the same origin, i.e., we need a 0 point of
% the ROI's projection that is invariant across iterations, so there is a
% unique correspondance of the grid points over time.
%cent = topo2inst(pointing0, gamma_topo(1), gamma_topo(2), target, sc, ...
%    inst, et);
%targetArea(:, 1) = targetArea(:, 1) - cent(1);
%targetArea(:, 2) = targetArea(:, 2) - cent(2);
% % Oriented area
% anglerot = -(angle - angle0);
% rotmat = [cosd(anglerot)   -sind(anglerot);
%           sind(anglerot)   cosd(anglerot)];
% % matrixGrid directions x and y
%[cxt, cyt] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
% orientedArea  = zeros(length(targetArea), 2);
% for j=1:length(targetArea)
%     orientedArea(j, :) = [cxt, cyt]' + rotmat*(targetArea(j, :)' - ...
%         [cxt, cyt]');
% end
% targetpshape = polyshape(orientedArea(:,1), orientedArea(:,2)); % polyshape

% Get grid shifting due to observation geometry update
updated_seed = topo2inst(gamma, cx, cy, target, sc, inst, et);
if ~isnan(updated_seed)
    shift = updated_seed - seed;
else
    shift = 0;
    updated_seed = seed;
    disp("Non existent gamma");
end
seed = updated_seed;
old_seed = old_seed + shift;

% Shift grid and tour
for i=1:size(inst_grid, 1)
    for j=1:size(inst_grid, 2)
        if ~isempty(inst_grid{i, j})
            inst_grid{i, j} = inst_grid{i, j} + shift;
        end
    end
end

for i=1:length(inst_tour)
    inst_tour{i} = inst_tour{i} + shift;
end

%% UPDATE GRID
% Update map by removing the previous element in the tour (next observation)
% Find which position does gamma occupy in this grid
ind_row = []; ind_col = [];
for i=1:numel(inst_grid)
    if ~isempty(inst_grid{i}) && norm(inst_grid{i} - old_seed) < 1e-3
        [ind_row, ind_col] = ind2sub(size(inst_grid), i);
        break;
    end
end
if ~isempty(old_seed), inst_grid{ind_row, ind_col} = []; end

% Enclose grid in a bigger matrix (with first and last rows and columns
% with NaN values, so we can explore neighbours adequately)
map = grid2map(inst_grid);

% Obtain the frontier tiles in the map: points that have less than 8
% neighbours in the grid
[frontier, indel] = getFrontierTiles(map);

% Update grid
openList = frontier; % open list starts as frontier set F
s = openList; % seeds: initial open list (frontier tiles)

while ~isempty(openList)

    % For each element in the openList array of grid points...
    o = openList{1}; % lon-lat coordinates of the observation point
    currind = indel{1}; % row-column indices of the observation point in
    % 'map'

    % Visited elements are deleted
    indel(1) = [];
    openList(1) = [];

    % Pre-allocate variables
    membershipChanged = false; % boolean variable that indicates if the
    % observation point has changed its membership in 'tour'
    insideTour = false; % boolean variable that indicates if the 
    % observation point belongs to 'tour'
    inSeed = false; % boolean variable that determines if the observation 
    % is in the seeding list

    % Previous check: is the current observation point inside the seeding
    % list? in case it is, let's analyze the neighbouring points (after
    % checking their membership in tour, next point)
    for i=1:length(s)
        if isequal(s{i}, o)
            inSeed = true;
            break;
        end
    end
    
    % Analyze the current element's membership in tour
    % Compute the current footprint's covered area
    fpshape = polyshape(o(1) + fpref.bvertices(:, 1), o(2) + fpref.bvertices(:, 2));
    inter = subtract(targetpshape, fpshape);
    areaI = area(inter);
    areaT = area(targetpshape);
    fpArea = area(fpshape);

    if (areaT - areaI)/fpArea >= epsilon % if the observation covers at
        % least a minimum roi's area...
        cin{end + 1} = o; % add it to the list of covering tiles
        cind{end +1} = currind;

        % Identify if the observation was already included in the
        % planned tour
        for i=1:length(inst_tour)
            if isequal(o, inst_tour{i})
                insideTour = true;
                break;
            end
        end
        if ~insideTour
            % If it wasn't included, then its membership changed
            membershipChanged = true;
        end
    else % otherwise (i.e., the observation's footprint falls outside the
        % roi)...
        cout{end + 1} = o; % add it to the list of disposal tiles

        % Identify if the observation was already included in the
        % planned tour
        for i=1:length(inst_tour)
            if isequal(o, inst_tour{i})
                insideTour = true;
                break;
            end
        end
        if insideTour
            % If it was included, then its membership changed
            membershipChanged = true;
        end
    end
    
    if inSeed || membershipChanged
        % Get the observation neighbor elements (diagonal and cardinal)
        [n, nind] = getNeighbours(o, currind, fpref.width, fpref.height, ...
            olapx, olapy, grid_dirx, grid_diry);

        % Check if the neighbors are already inside tour (in that case it
        % is not necessary to include them in the openlist for
        % re-evaluation)
        nindel = [];
        for i=1:length(n)
            inTour = false;
            for j=1:length(inst_tour)
                if norm(n{i} - inst_tour{j}) < 1e-5
                    inTour = true;
                    break;
                end
            end
            if inTour
                nindel = [nindel i];
            end
        end
        n(nindel) = [];
        nind(nindel) = [];

        for i=1:length(n)
            
            % Check if the neighbors are included in the cin or cout sets
            in1 = false; in2 = false;
            for j=1:length(cin)
                if norm(cin{j} - n{i}) < 1e-5
                    in1 = true;
                    break;
                end
            end
            for j=1:length(cout)
                if norm(cout{j} - n{i}) < 1e-5
                    in2 = true;
                    break;
                end
            end

            % If the neighbour node is not in the cin list nor in the
            % cout list... then add it to the openList for evaluation (if
            % not already included)
            if ~in1 && ~in2
                flag = false;
                for k=1:length(openList)
                    if norm(openList{k} - n{i}) < 1e-5
                        flag = true;
                        break;
                    end
                end
                if ~flag % if it's not in openList, then add it
                    openList{end+1} = n{i};
                    indel{end+1} = nind{i};
                end
            end
        end
    end
end

% Identify new tiles N = Cin - Tour
count = 0;
for i=1:length(cin)
    new = true; % let's assume that every element in 'cin' is new (and, 
    % thus, not included in 'tour' yet)
    for j=1:length(inst_tour)
        if norm(cin{i} - inst_tour{j}) < 1e-5
            new = false; % observation cin{i} is already included in 'tour'
            break;
        end
    end
    if new % if cin{i} is checked to be outside, include it in the new
        % tiles set
        count = count + 1;
        N{count} = cin{i};
        Nind{count} = cind{i};
    end
end

% Check that new identified tiles are not taboo (moving backwards in the
% coverage path)
ind_row = []; ind_col = [];
for i=2:numel(map)-1
    if ~any(isnan(map{i}), 'all') && norm(map{i} - seed) < 1e-3
        [ind_row, ind_col] = ind2sub(size(map), i);
        break;
    end
end
% Check that N is not coincident with old_seed...
for i=1:length(N)
    if norm(N{i} - old_seed) < 1e-5
        N(i) = []; Nind(i) = [];
        break;
    end
end
[N, Nind] = checkTaboo(N, Nind, map, ind_row, ind_col, sweepDir1, sweepDir2);

% Identify tiles to remove X = Cout - Tour
count = 0;
for i=1:length(cout)
    out = false; % let's assume that the disposable tiles are not 
    % included in 'tour'
    for j=1:length(inst_tour)
        if norm(cout{i} - inst_tour{j}) < 1e-5
            out = true; % observation cout{i} is included in 'tour'
            break;
        end
    end
    if out % if cout{i} is checked to be in 'tour', include it in the 
        % disposable tiles set
        count = count + 1;
        X{count} = cout{i};
    end
end

% Remove disposable tiles
map = removeTiles(map, X);

% Insert new tiles
map = insertTiles(map, N, Nind);

% Plot grid
% figure
% plot(targetpshape)
% hold on; axis equal;
% for i=1:size(grid, 1)
%     for j=1:size(grid, 2)
%         point = grid{i, j};
%         if ~isempty(point)
%             plot(point(1), point(2), 'b^')
%             hold on;
%         end
%     end
% end

% Boustrophedon decomposition
inst_grid = map2grid(map);
inst_tour = boustrophedon(inst_grid, sweepDir1, sweepDir2);

% Rotate back grid points
if ~isempty(inst_tour)
    topo_tour = inst2topo(inst_tour, cx, cy, target, sc, inst, et);
    % Remove empty elements from the tour, which may result from unobservable
    % regions within the planned plath
    emptyCells = cellfun(@isempty, topo_tour); % find indices of empty cells
    indEmpty = find(emptyCells);
    for k=1:length(indEmpty)
        emptyEl = inst_tour{indEmpty(k)};
        for i=1:size(map, 1)
            for j=1:size(map, 2)
                if ~isempty(map{i, j})
                    if vecnorm(map{i, j} - emptyEl) < 1e-5
                        map{i, j} = [];
                    end
                end
            end
        end
    end
    topo_tour(emptyCells) = []; % remove empty cells
    % Boustrophedon decomposition
    inst_grid = map2grid(map);
    inst_tour = boustrophedon(inst_grid, sweepDir1, sweepDir2);
    if ~isempty(inst_tour)
        seed = inst_tour{1};
    else
        seed = [];
    end
    %inst_tour(emptyCells) = []; % remove empty cells
    %seed = inst_tour{1};
    % for i=1:length(emptyCells)
    %     if ~emptyCells(i)
    %         seed = inst_tour{i};
    %         break;
    %     end
    % end
else
    seed = [];
    topo_tour = [];
end

end