function [n, fn] = getNeighboursInd(gamma, w, h, olapx, olapy, dx, dy, map, frontier)
% Given a point, this function outputs the 8 adjacent points.
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% Last Rev.:    06/2023
% 
% Usage:        n = getNeighbours(gamma, w, h, olapx, olapy, dx, dy)
%
% Inputs:
%   > gamma:        2D point to which the 8 adjacent points are needed
%   > w:            width (x-direction) of the grid spacing
%   > h:            height (y-direction) of the grid spacing
%   > olapx:        grid footprint overlap in the x direction,
%                   in percentage (width)
%   > olapy:        grid footprint overlap in the y direction,
%                   in percentage (height)
%   > dx:           vector that expresses the x-direction in the grid
%   > dy:           vector that expresses the y-direction in the grid
% 
% Outputs:
%   > n:            cell array with the 8 neighbouring points

% Pre-allocate variables
n = cell(8,1);
fn = cell(8,1);
ovlapx = olapx*w/100; ovlapy = olapy*h/100;

% Find gamma in the map
flag = false;
for i=1:size(map, 1)
    for j=1:size(map, 2)
        if norm(map{i, j} - gamma) < 1e-5
            indrow = i;
            indcol = j;
            flag = true;
            break;
        end
        if flag, break; end
    end
end

% Search neighbors of the given element in the map
el = map{indrow - 1, indcol + 1}; % northeast
if any(isnan(el)) || isempty(el)
    n{1} = gamma + (w-ovlapx)*dx +  (h-ovlapy)*dy; % northeast
elseif checkInFrontier(el, frontier)
    fn{1} = el;
end

el = map{indrow    , indcol + 1}; % east
if any(isnan(el)) || isempty(el)
    n{2} = gamma + (w-ovlapx)*dx ;                 % east
elseif checkInFrontier(el, frontier)
    fn{2} = el;
end

el = map{indrow + 1, indcol + 1}; % southeast
if any(isnan(el)) || isempty(el)
    n{3} = gamma + (w-ovlapx)*dx + (-h+ovlapy)*dy; % southeast
elseif checkInFrontier(el, frontier)
    fn{3} = el;
end

el = map{indrow - 1, indcol    }; % north
if any(isnan(el)) || isempty(el)
    n{4} = gamma + (h-ovlapy)*dy;                  % north
elseif checkInFrontier(el, frontier)
    fn{4} = el;
end

el = map{indrow + 1, indcol    }; % south
if any(isnan(el)) || isempty(el)
    n{5} = gamma + (-h+ovlapy)*dy;                 % south
elseif checkInFrontier(el, frontier)
    fn{5} = el;
end

el = map{indrow - 1, indcol - 1}; % northwest
if any(isnan(el)) || isempty(el)
    n{6} = gamma + (-w+ovlapx)*dx + (h-ovlapy)*dy; % northwest
elseif checkInFrontier(el, frontier)
    fn{6} = el;
end

el = map{indrow    , indcol - 1}; % west
if any(isnan(el)) || isempty(el)
    n{7} = gamma + (-w+ovlapx)*dx;                 % west
elseif checkInFrontier(el, frontier)
    fn{7} = el;
end

el = map{indrow + 1, indcol - 1}; % southwest
if any(isnan(el)) || isempty(el)
    n{8} = gamma + (-w+ovlapx)*dx +(-h+ovlapy)*dy; % southwest
elseif checkInFrontier(el, frontier)
    fn{8} = el;
end

end

function in = checkInFrontier(el, frontier)
in = false;
for i=1:length(frontier)
    if norm(el - frontier{i}) < 1e-5
        in = true;
        break;
    end
end
end