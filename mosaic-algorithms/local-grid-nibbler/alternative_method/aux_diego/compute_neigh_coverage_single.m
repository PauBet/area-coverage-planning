function neigh_coverage = compute_neigh_coverage_single(poly_roi,neighbours_polys)
     
    try  
       neigh_cov_poly = intersect(poly_roi,neighbours_polys); 
       footp_area = area(neighbours_polys);
       neigh_coverage = area(neigh_cov_poly)/footp_area;   
    catch
       neigh_coverage = 0;
    end    

end