function [N, X, map] = updateGridold(roi, tour, tour0, map, olapx, olapy, dirx, ...
    diry, fpref)
% Given a seed grid for the ROI discretization, this function analyzes its
% frontier points to search for new potential observations or disposable
% ones (no longer needed)
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% Last Rev.:    06/2023
% Pending revision
% 
% Usage:        [N, X, map] = updateGrid(roi, tour, map, olapx, ...
%                  olapy, dirx, diry, fpref)
%
% Inputs:
%   > roi:          matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D, in latitudinal 
%                   coordinates [º]
%       # roi(:,1) correspond to the x values of the vertices
%       # roi(:,2) correspond to the y values of the vertices         
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%   > map:          cell matrix of grid points. In order to avoid
%                   mapping boundaries, map is bounded by NaN rows and 
%                   columns (first and last)
%   > olapx:        grid footprint overlap in the x direction,
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction,
%                   in percentage (height)
%   > dirx:         vector that expresses the x-direction in the grid
%   > diry:         vector that expresses the y-direction in the grid
%   > fpref:        reference footprint that will be used to define the
%                   grid spacing of the roi's discretization. See footprint
%                   function for further details
% 
% Outputs:
%   > N:            cell matrix of the new observation points to be
%                   included in the tour
%   > X:            cell matrix of the disposable observation points to be
%                   removed from the tour
%   > map:          updated cell matrix of grid points (next observation
%                   point is removed from the map)

% Pre-allocate variables...
w = fpref.sizex; % width
h = fpref.sizey; % height
cin = {};  % set of new observations (inside or outside from 'tour')
cout = {}; % set of disposable observations (inside or outside from 'tour')
N = {};    % set of new tiles (outside from 'tour')
X = {};    % set of disposable tiles (inside of 'tour')

% Obtain the frontier tiles in the map: points that have less than 8
% neighbours in the grid
frontier = getFrontierTiles(map);

% Update grid
openList = frontier; % open list starts as frontier set F
s = openList; % seeds: initial open list
while ~isempty(openList)
    % For each element in the openList array of grid points...
    o = openList{1}; % lon-lat coordinates of the observation point
    openList(1) = []; % grid point checked
    membershipChanged = false; % boolean variable that indicates if the
    % observation point has changed its membership in 'tour'
    insideTour = false; % boolean variable that indicates if the 
    % observation point belongs to 'tour'
    inSeed = false; % boolean variable that determines if the observation 
    % is in the seeding list

    fpcurr = fptranslation(fpref, o(1), o(2)); % footprint allocation at 
    % current observation point

    % Previous check: is the current observation point inside the seeding
    % list? in case it is, let's move forward and analyze the neighbouring
    % points
    for i=1:length(s)
        if isequal(s{i}, o)
            inSeed = true;
            break;
        end
    end

    if ~inSeed % if the observation is not in the seeding list, then
        % analyze its membership in 'tour'

        % Compute the current footprint's covered area
        targetpshape = polyshape(roi(:,1), roi(:,2));
        fpshape = polyshape(fpcurr.bvertices(:, 1), fpcurr.bvertices(:, 2));
        inter = subtract(targetpshape, fpshape);
        areaI = area(inter);
        areaT = area(targetpshape);
        fpArea = area(fpshape);

        if (areaT - areaI)/fpArea >= 0.2 % if the observation covers at
            % least a minimum roi's area...
            cin{end + 1} = o; % add it to the list of covering tiles

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
    end
    
    if inSeed || membershipChanged
        % Get the observation neighbor elements (diagonal and cardinal)
        n = getNeighbours(o, w, h, olapx, olapy, dirx, diry);

        % Check if the neighbors are already inside tour (in that case it
        % is not necessary to include them in the openlist for
        % re-evaluation)
        indel = [];
        for i=1:length(n)
            inTour = false;
            for j=1:length(tour)
                if norm(n{i} - tour{j}) < 1e-5
                    inTour = true;
                    break;
                end
            end
            if inTour
                indel = [indel i];
            end
        end
        n(indel) = [];

        % Check if the neighbors are included in the cin or cout sets
        for i=1:length(n)
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
            % cout list... then add it to the openList for evaluation
            if ~in1 && ~in2
                flag = false;
                for k=1:length(openList)
                    if norm(openList{k} - n{i}) < 1e-5
                        flag = true;
                        break;
                    end
                end
                if ~flag, openList{end+1} = n{i}; end
            end
        end
    end
end

% Identify new tiles N = Cin - Tour
count = 0;
for i=1:length(cin)
    new = true; % let's assume that every element in 'cin' is new (and, 
    % thus, not included in 'tour' yet)
    for j=1:length(tour0)
        if norm(cin{i} - tour0{j}) < 1e-5
            new = false; % observation cin{i} is already included in 'tour'
            break;
        end
    end
    if new % if cin{i} is checked to be outside, include it in the new
        % tiles set
        count = count + 1;
        N{count} = cin{i};
    end
end

% Identify tiles to remove X = Cout - Tour
count = 0;
for i=1:length(cout)
    out = false; % let's assume that the disposable tiles are not included
    % in 'tour'
    for j=1:length(tour0)
        if norm(cout{i} - tour0{j}) < 1e-5
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

% Update map by removing the first element in the tour (next observation)
flag = false;
for i=2:size(map, 1)
    for j=2:size(map, 2)
        if norm(map{i, j} - tour{1}) < 1e-5
            map{i, j} = [NaN, NaN];
            flag = true;
            break;
        end
    end
    if flag, break; end
end


% figure
% hold on; box on;
% for i=1:size(map, 1)
%     for j=1:size(map, 2)
%         plot(map{i, j}(1), map{i, j}(2), 'b*')
%     end
% end
end