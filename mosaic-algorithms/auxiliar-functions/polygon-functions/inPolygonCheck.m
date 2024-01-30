function in = inPolygonCheck(xi, yi, x, y, xq, yq)

%in = inpolygon(xi, yi, x, y);
in = inpolygon(xi, yi, xq, yq);
if ~in
    in = inpolygon(x, y, xq, yq);
    for k=1:length(in)
        if in(k)
            in = 1;
            break;
        end
    end
end

end