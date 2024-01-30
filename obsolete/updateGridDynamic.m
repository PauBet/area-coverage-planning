function [N, X, map] = updateGridDynamic(roi, tour, map, ...
    t, sc, inst, target, theta, olapx, olapy, dirx, diry, fprinti)

%frontier = getFrontierTiles(tour, map);
frontier = tour{1};
bbox  = smallestBoundingBox(fprinti.bvertices(:, 1), fprinti.bvertices(:, 2));
w = bbox.size1;
h = bbox.size2;

polygonArea = polyshape(roi(:,1), roi(:,2));

%% temporary figure
figure
plot(polygonArea)
hold on;
for i=1:size(map, 1)
    for j=1:size(map, 2)
        if ~any(isnan(map{i, j})) && ~any(isempty(map{i, j}))
            plot(map{i, j}(1), map{i, j}(2), 'r*')
        end
    end
end
for i=2:length(tour)
    plot([tour{i-1}(1) tour{i}(1)], [tour{i-1}(2) tour{i}(2)],...
            'w-', 'linewidth', 1)
end
drawnow

cin = {}; % closed list
cout = {};% closed list
openList{1} = frontier'; % open list starts as frontier set F
if isempty(openList)
    %[openList(:,1), openList(:,2)] = discretize(bbox, sizex, sizey, overlapx, overlapy);
end
s = openList; % seeds: initial open list
while ~isempty(openList)
    o = openList{1}';
    openList(1) = [];
    membershipChanged = false; % boolean variable that indicates if the
    % observation point has changed its membership to 'tour'
    insideTour = false; % boolean variable that indicates if the 
    % observation point belongs to 'tour'

    % Observation footprint at time t
    fprintobs = footprint(o(1), o(2), t, inst, sc, target, theta);

    % Compute the current footprint's covered area
    targetpshape = polyshape(roi(:,1), roi(:,2));
    fpshape = polyshape(fprintobs.bvertices(:, 1), fprintobs.bvertices(:, 2));
    inter = subtract(targetpshape, fpshape);
    areaI = area(inter);
    areaT = area(targetpshape);
    fpArea = area(fpshape);
    
    if (areaT - areaI)/fpArea >= 0.1 % if the observation covers at least
        % a minimum roi's area...
        cin{end + 1} = o; % add it to the list of covering tiles

        % Identify if the observation was already included in the planned
        % tour
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

        % Identify if the observation was already included in the planned
        % tour
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

    % Determine if the observation is in the seeding list
    inSeed = false;
    if ~membershipChanged % if membership already changed, then there is no
        % need to check this...
        for i=1:length(s)
            if isequal(s{i}', o)
                inSeed = true;
                break;
            end
        end
    end

    if membershipChanged || inSeed
        % Get the observation neighbor elements (diagonal and cardinal)
        n = getNeighbors(o', w, h, olapx, olapy, dirx, diry);

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
            % cout list... then add it to the openList
            if ~in1 && ~in2
                flag = false;
                for k=1:length(openList)
                    if norm(openList{k} - n{i}') < 1e-5
                        flag = true;
                        break;
                    end
                end
                if ~flag, openList{end+1} = n{i}'; end
            end
        end
    end
end

% Identify new tiles N = Cin - Tour
count = 0;
N = {};
for i=1:length(cin)
    new = true;
    for j=1:length(tour)
        if norm(cin{i} - tour{j}) < 1e-5
            new = false;
            break;
        end
    end
    if new
        count = count + 1;
        N{count} = cin{i};
        new = true;
    end
end

% Identify tiles to remove
count = 0;
out = false;
X = {};
for i=1:length(cout)
    for j=1:length(tour)
        if norm(cout{i} - tour{j}) < 1e-5
            out = true;
            break;
        end
    end
    if out
        count = count + 1;
        X{count} = cout{i};
        out = false;
    end
end

flag = false;
for i=2:size(map, 1)
    for j=2:size(map, 2)
        if norm(map{i, j} - tour{1}') < 1e-5
            map{i, j} = [NaN, NaN];
            flag = true;
            break;
        end
    end
    if flag, break; end
end

%% temporary figure
% figure
% plot(polygonArea)
% hold on;
% for i=1:length(cin)
%     plot(cin{i}(1), cin{i}(2), 'g*')
% end
% for i=1:length(cout)
%     plot(cout{i}(1), cout{i}(2), 'k*')
% end
% for i=1:length(N)
%     plot(N{i}(1), N{i}(2), 'r*')
% end
% for i=1:length(X)
%     plot(X{i}(1), X{i}(2), 'b*')
% end

end