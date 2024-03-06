function [A, fplist, time, rotmatrix, poly_roi, s_area] = target_coverage(target, time, inst, sc, target_body, poly_roi, A, fplist)

    [target_fp, rotmatrix]= footprint(time, inst, sc, target_body, 'lowres', target(1), target(2), 1);
                
    poly_target_footprint = polyshape(target_fp.bvertices);

    % Compute the coverage of the current footprint over the real ROI    
    real_coverage = intersect(poly_roi,poly_target_footprint);
    fpcoverage = area(real_coverage)/area(poly_target_footprint);
    % Compute the remaining area after substracting the footprint to the
    % ROI
    poly_roi = subtract(poly_roi,poly_target_footprint);
    s_area = area(poly_roi);
                
    % If the actual footprint covers part of the real ROI save it and move
    % to the next instant
    
    if real_coverage.NumRegions ~= 0 
        if fpcoverage >= 0.2
            A{end+1} = [target(1), target(2)];
            fplist(end + 1) = target_fp;
        end
    end

end