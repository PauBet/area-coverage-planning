function [next_target, target_fp, poly_target_footprint, neigh_indexes, rotmatrix, fpcoverage, count] = select_next_target(neigh_indexes,neighbours,poly_roi,ptime,footprint_func, rot_func, count)

% From the valid neighbours this function selects the next target to be evaluated. 
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date:        08/2023
% Version:      1
% Last update: 03/2024
%
% Usage:  [next_target, target_fp, poly_target_footprint, neigh_indexes, rotmatrix, fpcoverage, count] = select_next_target(neigh_indexes,neighbours,poly_roi,...
%                                                                                                             ptime,footprint_func, rot_func, count)
%
% Inputs: 
%    > neigh_indexes:           List of integers with the indices of the
%                               valid neighbours. 
%    > neighbours:              Coordinates of the valid neighbous
%    > poly_roi:                Polyshape object containing the remaining
%                               ROI.
%    > ptime:                   List of possible imaging times of the
%                               neighbours.
%    > footprint_func           Function handle containing the footprint
%                               function, so that only time and target
%                               coodinates are needed.
%    > rot_func                 Function handle containing the
%                               rotation/pointing function, so that only time and target
%                               coodinates are needed.
%    > count:                   Counter of fooptrint_func calls (analytical purposes)
% Outputs:
%    > next_target:             Index of the neighbour chosen as the next target
%    > target_fp:               Footprint struct of the selected neighbour
%    > poly_target_footprint: 
%    > neigh_indexes:
%    > rotmatrix:
%    > fpcoverage:
%    > count                    Counter of fooptrint_func calls (analytical purposes)            

    % Preallocate in case no valid target is found
    next_target = 0;
    target_fp = [];
    rotmatrix = [];
    poly_target_footprint = [];
    fpcoverage = [];
    
    % Define a minimum acceptable coverage area (less that the value
    % required for adding it to the final list).
    min_area = 0.05;

    n = 1;

    for i = neigh_indexes
        
        % Following the order of neigh_indexes, compute for each valid
        % neighbour its fooptrint coverage, and check if it meets the
        % minimum criteria.

        [neigh_poly{i}, neigh_rotation_mat(:,:,i), neigh_fp{i}] =  get_poly_rot(neighbours(i,:), ptime(i), footprint_func, rot_func);
        neigh_coverage(i) = get_fp_coverage(poly_roi,neigh_poly{i});
        
        count = count + 1;
        
        % If the neighbour meets the minimum assign its index and related
        % variables to the ones of the next target

        if neigh_coverage(i)>min_area
            next_target = i;
            if n>1
                % If n>1 a edge of the ROI has been reached, and the path
                % will start going in the opposite direction. Therefore the
                % list of indexes has to be updated. 
                opposite = get_opposite_neigh(neigh_indexes(1));

                % Now the opposite of the previous main direction will be
                % checked, and the other neighbours will be checked in
                % reverted order (if before (east, south east, south, south west) now 
                %                 (west, south west, south, south east)).
                neigh_indexes = [opposite, flip(neigh_indexes(2:end))];
            end 
            
            % Assign outputs
            target_fp = neigh_fp{i};
            rotmatrix = neigh_rotation_mat(:,:,i);
            poly_target_footprint = neigh_poly{i};
            fpcoverage = neigh_coverage(i);

            break
        end    
        n = n + 1;
    end

end