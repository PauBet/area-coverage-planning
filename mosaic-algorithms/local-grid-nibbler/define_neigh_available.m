function neigh_indexes = define_neigh_available(dist_ind, neigh_coverage)
  
% From the coverage of the 8 possible neighbours of the current target, and knowing 
% which neighbour aligns better with the precomputed coverage path, it
% identifies a list of 4 neighbours that will be used to cover the ROI, the
% other 4 are located at unvalid locations that are not contained in the
% coverage path, therefore, they are discarded to reduce computational
% cost.
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date:        08/2023
% Version:      1
% Last update: 03/2024
%
% Usage:       neigh_indexes = define_neigh_available(dist_ind, neigh_coverage)
%
% Inputs: 
%    > dist_ind:       Index of the neighbour that aligns with the coverage
%                      path at the current point.
%    > neigh_coverage: List with the coverage of the 8 neighbours
% Outputs:
%    > neigh_indexes: List with the indexes of the valid neighbours


            % Neigh coverage sort           
            [~, nei_ind] = sort(neigh_coverage, 'descend');


            % Select which neighbours will be used in the following
            % iterations

            if nei_ind(1) ~= dist_ind(1)
                max_cn = nei_ind(1);
            else
                max_cn = nei_ind(2);
            end

            nei_grid = [8,1,2; 2,3,4; 6,5,4; 6,7,8];
            [possible_lines_index, ~] = find(nei_grid == max_cn);
            possible_lines = nei_grid(possible_lines_index,:);
            
            % Discard line with next target
            [invalid_lines_index, ~] = find(possible_lines == dist_ind(1));
            possible_lines(invalid_lines_index,:) = [];

            if size(possible_lines,1) > 1
                cov_sum(1) = sum(neigh_coverage(possible_lines(1,:)));
                cov_sum(2) = sum(neigh_coverage(possible_lines(2,:)));
                if cov_sum(1)>cov_sum(2)
                    possible_lines(2,:) = [];
                else
                    possible_lines(1,:) = [];
                end    
            end
            
            neigh_indexes = [dist_ind(1), possible_lines]; % Neighbours to be visited in next iterations
            if abs(neigh_indexes(1)-neigh_indexes(2)) ~= 1
                neigh_indexes = [dist_ind(1), flip(possible_lines)];
            end    
end