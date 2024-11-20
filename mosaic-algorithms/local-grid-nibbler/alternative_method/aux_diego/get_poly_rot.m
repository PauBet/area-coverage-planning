function [footprint_poly, rotation_mat, target_fp] =  get_poly_rot(target, time, footprint_func)
    
    %For a given target in lon/lat provides a footprint in polyshape format
    %and struct format
    %and a rotation matrix needed for later
    [target_fp, rotation_mat] = footprint_func(target(1,1), target(1,2),time);
    target_footprint = target_fp.bvertices;
    try
        footprint_poly = polyshape(target_footprint);
    catch
        footprint_poly = 0;
    end    

end