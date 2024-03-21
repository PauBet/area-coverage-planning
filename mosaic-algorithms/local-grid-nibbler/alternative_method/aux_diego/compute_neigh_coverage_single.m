function fp_coverage = get_fp_coverage(poly_roi,fp_poly)
     
    try  
       intersec = intersect(poly_roi,fp_poly); 
       intersec_area = area(intersec);
       fp_coverage = area(fp_poly)/intersec_area;   
    catch
       fp_coverage = 0;
    end    

end