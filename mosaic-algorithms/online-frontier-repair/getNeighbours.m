function [n, nind] = getNeighbours(gamma, ind, w, h, olapx, olapy, dx, dy)
% Given a point, this function outputs the 8 adjacent neighbours and their
% index location in the grid
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% Last Rev.:    06/2023
% 
% Usage:        [n, nind] = getNeighbours(gamma, ind, w, h, olapx, ...
%                   olapy, dx, dy)
%
% Inputs:
%   > gamma:        2D point to which the 8 adjacent points are needed
%   > ind:          index position of gamma in the grid
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
%   > nind:         cell array with the index position of the neighbours
%                   in the grid

% Pre-allocate variables
n = cell(8,1);
ovlapx = olapx*w/100; ovlapy = olapy*h/100;

% Cardinal and diagonal neighbouring points
n{1} = gamma + (w-ovlapx)*dx +  (h-ovlapy)*dy; % northeast
n{2} = gamma + (w-ovlapx)*dx ;                 % east
n{3} = gamma + (w-ovlapx)*dx + (-h+ovlapy)*dy; % southeast
n{4} = gamma + (h-ovlapy)*dy;                  % north
n{5} = gamma + (-h+ovlapy)*dy;                 % south
n{6} = gamma + (-w+ovlapx)*dx + (h-ovlapy)*dy; % northwest
n{7} = gamma + (-w+ovlapx)*dx;                 % west
n{8} = gamma + (-w+ovlapx)*dx +(-h+ovlapy)*dy; % southwest

% Cardinal and diagonal neighbouring points grid indices
nind{1} = ind + [-1  1];
nind{2} = ind + [ 0  1];
nind{3} = ind + [ 1  1];
nind{4} = ind + [-1  0];
nind{5} = ind + [ 1  0];
nind{6} = ind + [-1 -1];
nind{7} = ind + [ 0 -1];
nind{8} = ind + [ 1 -1];

end