function bbox = computeBoundingBox(vertices)
% [Description]

%%
maxlon = -inf;
maxlat = -inf;
minlat = inf;
minlon = inf;

for i=1:length(vertices)
    if vertices(i,1) > maxlon
        maxlon = vertices(i,1);
    end
    if vertices(i,1) < minlon
        minlon = vertices(i,1);
    end

    if vertices(i,2) > maxlat
        maxlat = vertices(i,2);
    end
    if vertices(i,2) < minlat
        minlat = vertices(i,2);
    end
end

%%
bbox.maxlon = maxlon;
bbox.minlon = minlon;
bbox.maxlat = maxlat;
bbox.minlat = minlat;
end