function [N, X, F] = updateGrid(vertices, tour, F, map)
polygonArea = polyshape(vertices(:,1), vertices(:,2));
cin = {}; % closed list
cout = {};% closed list
openList = F; % open list starts as frontier set F
if isempty(openList)
    %[openList(:,1), openList(:,2)] = discretize(bbox, sizex, sizey, overlapx, overlapy);
end
s = openList; % seeds: initial open list
openList = cell2mat(openList');
while ~isempty(openList)
    o = openList(1,:);
    openList(1,:) = [];
    membershipChanged = false;
    inside = false;
    if inpolygon(o(1), o(2), vertices(:,1), vertices(:,2))
        cin{end + 1} = o;
        for i=1:length(tour)
            if isequal(o, tour{i})
                inside = true;
                break;
            end
        end
        if ~inside
            %tour{i} = [];
            membershipChanged = true;
        end
    else
        cout{end + 1} = o;
        for i=1:length(tour)
            if isequal(o, tour{i})
                inside = true;
                break;
            end
        end
        if inside
            %tour{end + 1} = o;
            membershipChanged = true;
        end
    end
    in = false;
    for i=1:length(s)
        if isequal(s{i}, o)
            in = true;
            break;
        end
    end

    if membershipChanged || in
        for i=1:size(map,1)
            for j=1:size(map,2)
                if isequal(map{i,j}, o)
                    indrow = i;
                    indcol = j;
                    break;
                end
            end
        end
        n = getNeighbors(indrow, indcol, map);
        in1 = false; in2 = false;
        for i=1:length(n)
            for j=1:length(cin)
                if isequal(cin{j}, n{i})
                    in1 = true;
                    break;
                end
            end
            for j=1:length(cout)
                if isequal(cout{j}, n{i})
                    in2 = true;
                    break;
                end
            end
            %if xor(~in1, ~in2)
            if ~in1 && ~in2
                openList(end+1, :) = cell2mat(n(i));
            end
        end
    end
end

% Identify new tiles N = Cin - Tour
count = 0;
new = true;
N = {};
for i=1:length(cin)
    for j=1:length(tour)
        if isequal(cin{i}, tour{j})
            new = false;
            break;
        end
    end
    if new
        count = count + 1;
        N{count} = cin{i};
    end
end

% Identify tiles to remove
count = 0;
out = false;
X = {};
for i=1:length(cout)
    for j=1:length(tour)
        if isequal(cout{i}, tour{j})
            out = true;
            break;
        end
    end
    if out
        count = count + 1;
        X{count} = cout{i};
    end
end

%N = cin - tour; 
F = {};

for i=1:length(tour) % is it the updated tour or the previous one?
    tile = tour{i};
    for j=1:size(map,1)
        exit = false;
        for k=1:size(map,2)
            if isequal(tile, map{j,k})
                exit = true;
                break;
            end
        end
        if exit
            break;
        end
    end
    n = getNeighbors(j, k, map);
    aux_tour = zeros(length(tour),2);
    aux_n = zeros(length(n),2);
    for t=1:length(tour)
        aux_tour(t,:) = tour{t};
    end
    for nn=1:length(n)
        aux_n(nn,:) = n{nn};
    end
    if length(intersect(aux_n, aux_tour, 'rows')) < 8
        F{end+1} = tile;
    end
end

end