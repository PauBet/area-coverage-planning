function ax = fp2mapplot(ax, target, fplist, varargin)

% Pre-allocate variables
filename = strcat(lower(target),'-map.jpg');
if nargin > 3 
    facecolor = varargin{1};
    edgecolor = varargin{2};
else
    facecolor = 'c';
    edgecolor = 'c';
end

% Figure
if isempty(ax), ax = mapPlot(filename); end
title('Coverage map')
set(ax, 'FontSize', 20)
% Plot footprint list
for i=1:length(fplist)
    if ~isempty(fplist(i).bvertices)
        poly = polyshape(fplist(i).bvertices(:, 1), fplist(i).bvertices(:, 2));
        plot(poly, 'FaceColor', facecolor, 'FaceAlpha', 0.2, 'EdgeColor', edgecolor);
    end
end
end