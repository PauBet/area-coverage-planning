function opposite = get_opposite_neigh(neigh)
% Since the neighbours are included in a 3x3 grid each neighbour that in
% not in a diagonal position, has a opposite neighbour. When the algorithms
% reaches an edge of the ROI, it has to move down or up, and start moving
% in the opposite direction, which is given by this function.
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date:        08/2023
% Version:      1
% Last update: 03/2024
%
% Usage:       opposite = get_opposite_neigh(neigh)
%
% Inputs: 
%    > neigh: Index of the neighbour whose opposite we are looking for.
% Outputs:
%    > opposite: Index of the opposite neighbour.

    if neigh == 1
        opposite = 5;
    elseif neigh == 3
        opposite = 7;
    elseif neigh == 5
        opposite = 1;
    elseif neigh == 7
        opposite = 3;
    end     

end
