function best = checkNeighbors(gamma, map, center)
% [Description]

%%
for i=1:size(map,1)
    for j=1:size(map,2)
        if norm(map{i,j} - gamma) < 1e-5
            indrow = i;
            indcol = j;
            break;
        end
    end
end
c = getNeighbors(indrow, indcol, map, 'cardinal');
d = getNeighbors(indrow, indcol, map, 'diagonal');
%% argmax
bestScore = 0;
for i=1:length(c)
    currEl = c{i};
    currScore = abs(score(currEl, center));
    if currScore >= bestScore
        bestScore = currScore;
        best = currEl;
    end
end
%%
if 1
    x = 1.1; % bias against diagonal neighbors
    bestScore = bestScore*x;
    for i=1:length(d)
        currScore = score(d{i}, center);
        if 1
             best = d{i};
        elseif currScore > bestScore
            bestScore = currScore;
            best = d{i};
        end
    end
end

end