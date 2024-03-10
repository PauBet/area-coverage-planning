function [gamma, grid, tour, tour_topo] = updateGrid2(roi, ...
    tour, grid, grid_dirx, grid_diry, angle0, cx, cy, olapx, olapy, angle, ...
    indir1, indir2, gamma, gamma0, gamma_topo, et, inst, sc, target)
%(roi, tour, ...
%    tour_topo, olapx, olapy, indir1, indir2, grid_topo, origin_topo, ...
%    old_origin_topo, grid, gamma, gamma0, et, inst, sc, target, angle)

% Pre-allocate variables...
persistent fpref;
persistent pointing0;
persistent dir1;
persistent dir2;
if isempty(dir1)
    dir1 = indir1; dir2 = indir2;
end
if isempty(pointing0), [pointing0(1), pointing0(2)] = deal(cx, cy); end
cin = {}; % set of new observations (inside or outside from 'tour')
cind= {}; % map indices of the new observation (potential placement in the
% map)
cout = {}; % set of disposable observations (inside or outside from 'tour')
N    = {}; % set of new tiles (outside from 'tour')
Nind = {}; % map indices of the new tiles
X    = {}; % set of disposable tiles (inside of 'tour')
epsilon = 0.2;

% Build reference tile (it's always going to be the same in the subsequent
% calls)
if isempty(fpref)
    % We need to craft the tile reference that we are going to use 
    % throughout the heruristics operations. To avoid repetitions, we
    % declare it as a persistent variable that is only going to be
    % calculated in the first call
    [~, ~, ~, bounds] = ...
        cspice_getfov(cspice_bodn2c(inst), 4); % get fovbounds in the instrument's reference frame
    [~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
    fpref.width = width;
    fpref.height = height;
    fpref.angle = angle0;
    xlimit = [-width/2 width/2]; ylimit = [-height/2 height/2];
    xbox = xlimit([1 1 2 2 1]); ybox = ylimit([1 2 2 1 1]);
    fpref.bvertices(:, 1) = xbox;
    fpref.bvertices(:, 2) = ybox;
end

% Project ROI topographical coordinates to instrument's focal plane
targetArea = topo2inst(roi, cx, cy, target, sc, inst, et);
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
targetpshape = polyshape(targetArea(:,1), targetArea(:,2)); % polyshape
% Oriented area
% anglerot = -(angle - angle0);
% rotmat = [cosd(anglerot)   -sind(anglerot);
%           sind(anglerot)   cosd(anglerot)];
% % matrixGrid directions x and y
% [cxt, cyt] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
% orientedArea  = zeros(length(targetArea), 2);
% for j=1:length(targetArea)
%     orientedArea(j, :) = [cxt, cyt]' + rotmat*(targetArea(j, :)' - ...
%         [cxt, cyt]');
% end
% targetpshape = polyshape(orientedArea(:,1), orientedArea(:,2)); % polyshape

% Get grid shifting due to observation geometry update
updated_gamma = topo2inst(gamma_topo, cx, cy, target, sc, inst, et);
shift = updated_gamma - gamma;
gamma = updated_gamma;
gamma0 = gamma0 + shift;

% Shift grid, tour, gamma and gamma0
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        if ~isempty(grid{i, j})
            grid{i, j} = grid{i, j} + shift;
        end
    end
end

for i=1:length(tour)
    tour{i} = tour{i} + shift;
end

%% UPDATE GRID
% Update map by removing the previous element in the tour (next observation)
% Find which position does gamma occupy in this grid
ind_row = []; ind_col = [];
for i=1:numel(grid)
    if ~isempty(grid{i}) && norm(grid{i} - gamma0) < 1e-3
        [ind_row, ind_col] = ind2sub(size(grid), i);
        break;
    end
end
if ~isempty(gamma0), grid{ind_row, ind_col} = []; end

% Enclose grid in a bigger matrix (with first and last rows and columns
% with NaN values, so we can explore neighbours adequately)
map = grid2map(grid);

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
    
    % Tile allocation at the current observation point
    %fpcurr = fptranslation(fpref, o(1), o(2));

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
        for i=1:length(tour)
            if isequal(o, tour{i})
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
        for i=1:length(tour)
            if isequal(o, tour{i})
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
            for j=1:length(tour)
                if norm(n{i} - tour{j}) < 1e-5
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
    for j=1:length(tour)
        if norm(cin{i} - tour{j}) < 1e-5
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
    if ~any(isnan(map{i}), 'all') && norm(map{i} - gamma) < 1e-3
        [ind_row, ind_col] = ind2sub(size(map), i);
        break;
    end
end
[N, Nind] = checkTaboo2(N, Nind, map, ind_row, ind_col, dir1, dir2);

% Identify tiles to remove X = Cout - Tour
count = 0;
for i=1:length(cout)
    out = false; % let's assume that the disposable tiles are not 
    % included in 'tour'
    for j=1:length(tour)
        if norm(cout{i} - tour{j}) < 1e-5
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
[~, map] = removeTiles(tour, map, X);

% Insert new tiles
[~, map] = insertTiles(tour, map, N, Nind);

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
grid = map2grid(map);
tour = boustrophedon(grid, dir1, dir2);

% Rotate back grid points
%grid_rot = rotGrid([cxt cyt], grid, angle, angle0);
%tour_rot = rotGrid([cxt cyt], tour, angle, angle0);
if ~isempty(tour)
    gamma = tour{1};
    % Transform coordinates
    % for i=1:size(grid_rot, 1)
    %     for j=1:size(grid_rot, 2)
    %         if ~isempty(grid_rot{i, j}) 
    %             grid_rot{i, j} = grid_rot{i, j} + cent;
    %         end
    %     end
    % end
    % for i=1:length(tour_rot)
    %         if ~isempty(tour_rot{i})
    %             tour_rot{i} = tour_rot{i} + cent;
    %         end
    % end
    %grid_topo = inst2topo(grid_rot, gamma_topo(1), gamma_topo(2), target, sc, inst, et);
    %tour_topo = inst2topo(tour_rot, gamma_topo(1), gamma_topo(2), target, sc, inst, et);
    grid_topo = inst2topo(grid, cx, cy, target, sc, inst, et);
    tour_topo = inst2topo(tour, cx, cy, target, sc, inst, et);
    indel = [];
    for i=1:length(tour_topo)
        if isempty(tour_topo{i}), indel = [indel i]; end
    end
    tour_topo(indel) = [];
    %grid_topo = rotGrid(pointing0, grid_topo, angle, angle0);
    %tour_topo = rotGrid(pointing0, tour_topo, angle, angle0);
else
    gamma = [];
    grid_topo = [];
    tour_topo = [];
end

end

function matrixGrid = rotGrid(cent, grid, inangle, angle0)
    angle = inangle - angle0;
    angle = -angle;
    matrixGrid = grid;
    if angle ~= 0
        rotmat = [cosd(angle)   -sind(angle);
            sind(angle)   cosd(angle)];
        cx = cent(1); cy = cent(2);

        % Sort and rotate the grid points and insert them in the grid matrix
        matrixGrid = cell(size(grid));
        for i=1:size(grid, 1)
            % We will sweep across the grid by, first, latitude and, second,
            % longitude
            for j=1:size(grid, 2)
                if ~isempty(grid{i, j})
                    lon = grid{i, j}(1); lat = grid{i, j}(2);
                    matrixGrid{i, j} = ([cx, cy]' + transpose(rotmat)*...
                        ([lon lat] - [cx, cy])')'; % rotate values to the original roi
                    % orientation
                end
            end
        end
    end
end