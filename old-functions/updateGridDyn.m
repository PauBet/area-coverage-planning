function [N, X, map] = updateGridDyn(roi, tour, map, dx, dy, olapx, olapy, dirx, ...
    diry, fpref)

% Pre-allocate variables...
w = dx; % width
h = dy; % height
cin = {};  % set of new observations (inside or outside from 'tour')
cout = {}; % set of disposable observations (inside or outside from 'tour')
N = {};    % set of new tiles (outside from 'tour')
X = {};    % set of disposable tiles (inside of 'tour')
ovlapx = olapx*w/100; ovlapy = olapy*h/100;

% Analyze inner points of the map -> these points shall not be modified
% Obtain the frontier tiles in the map: points that have less than 8
% neighbours in the grid
[frontier, inner] = getFrontierTiles(map, tour); % frontier and inner
% points are sorted in order of appearance in the 'tour'
openList = frontier; % open list starts as frontier set F
s = openList; % seeds: initial open list
visitedFrontier = {};

% Sweep across the planned tour
while ~isempty(openList)
    % For each observation in the planned tour
    o = openList{1}; % lon-lat coordinates of the observation point
    openList(1) = []; % grid point checked
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

    if ~inSeed
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
        else % otherwise (i.e., the observation's footprint falls outside the
            % roi)...
            cout{end + 1} = o; % add it to the list of disposal tiles
        end
    end

    % Temporary figure
    % figure
    % hold on;
    % for i=1:length(tour)
    %     plot(tour{i}(1), tour{i}(2), 'bo', 'MarkerSize', 15)
    % end
    % plot(o(1), o(2), 'k*')

    if inSeed

        % Get the observation neighbor elements (diagonal and cardinal)
        [n, fn] = getNeighboursInd(o, fpref.sizex, fpref.sizey, olapx, olapy, dirx, diry, map, frontier);
        % Note: there is no need to check if 'n' belong to 'tour' because
        % the function already outputs elements outside from it (inner
        % points are not considered as neighbours and frontier points are
        % saved in fn)

        for i=1:length(n)

            % Case 1: For each outside neighbour, this is, a neighbouring
            % element of the grid that is not included in the tour,
            % consider it for inclusion (openList)
            if ~isempty(n{i})
                % Check if neighbours should be included in the openList
                % (no previously visited neighbours nor duplicates in the 
                % set)

                in1 = false; in2 = false;
                for j=1:length(cin)
                    if norm(cin{j} - n{i}) < 1e-5
                        in1 = true;
                        break;
                    end
                end
                for j=1:length(cin)
                    if norm(cout{j} - n{i}) < 1e-5
                        in1 = true;
                        break;
                    end
                end

                % If the neighbour node is not in the cin list nor in the
                % cout list... then add it to the openList for evaluation
                if ~in1 && ~in2
                    openList{end+1} = n{i};
                end
            
            % Case 2: For each frontier tile, re-evaluate with the current
            % grid spacing (fpref)
            elseif ~isempty(fn{i})

                % Check that the frontier tile has not been previously
                % visited
                flag = false;
                for k=1:length(visitedFrontier)
                    if norm(visitedFrontier{k} - fn{i}) < 1e-5
                        flag = true;
                        break;
                    end
                end
                if flag, continue; end

                % Find the frontier element position in the tour and the
                % map
                flag = false;
                for ii=2:size(map, 1)
                    for jj=2:size(map, 2)
                        if norm(map{ii, jj} - fn{i}) < 1e-5
                            flag = true;
                            break;
                        end
                    end
                    if flag, break; end
                end

                for tt=1:length(tour)
                    if norm(tour{tt} - fn{i}) < 1e-5
                        break;
                    end
                end

                for oo=1:length(openList)
                    if norm(openList{oo} - fn{i}) < 1e-5
                        break;
                    end
                end

                for ss=1:length(s)
                    if norm(s{ss} - fn{i}) < 1e-5
                        break;
                    end
                end

                for ff=1:length(frontier)
                    if norm(frontier{ff} - fn{i}) < 1e-5
                        break;
                    end
                end

                neww = fpref.sizex;
                newh = fpref.sizey;
                % Correct neighbour position
                switch i
                    case 1
                        fn{1} = o + (neww-ovlapx)*dirx +  (newh-ovlapy)*diry; % northeast
                    case 2
                        fn{2} = o + (neww-ovlapx)*dirx ;                 % east
                    case 3
                        fn{3} = o + (neww-ovlapx)*dirx + (-newh+ovlapy)*diry; % southeast
                    case 4
                        fn{4} = o + (newh-ovlapy)*diry;                  % north
                    case 5
                        fn{5} = o + (-newh+ovlapy)*diry;                 % south
                    case 6
                        fn{6} = o + (-neww+ovlapx)*dirx + (newh-ovlapy)*diry; % northwest
                    case 7
                        fn{7} = o + (-neww+ovlapx)*dirx;                 % west
                    case 8
                        fn{8} = o + (-neww+ovlapx)*dirx +(-newh+ovlapy)*diry; % southwest
                end

                map{ii, jj} = fn{i};
                tour{tt} = fn{i};
                openList{oo} = fn{i};
                s{ss} = fn{i};
                frontier{ff} = fn{i};

                visitedFrontier{end+1} = fn{i};

            end
        end

    end
    
    % for i=1:length(tour)
    %     plot(tour{i}(1), tour{i}(2), 'cp', 'MarkerSize', 10)
    % end
    % for i=1:length(n)
    %     if ~isempty(n{i})
    %         plot(n{i}(1), n{i}(2), 'r^', 'MarkerSize', 15)
    %     end
    %     if ~isempty(fn{i})
    %         plot(fn{i}(1), fn{i}(2), 'g^', 'MarkerSize', 15)
    %     end
    % end

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
end