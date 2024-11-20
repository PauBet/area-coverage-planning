function neigh_indexes = define_neigh_available(dist_ind, neigh_coverage)
            
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