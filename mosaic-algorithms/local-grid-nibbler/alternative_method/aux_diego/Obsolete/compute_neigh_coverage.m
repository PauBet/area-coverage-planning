function neigh_coverage = compute_neigh_coverage(poly_roi,neighbours_polys,neigh_indexes)
    
    neigh_coverage = zeros(1,8);

    for i = neigh_indexes       
        try  
            neigh_cov_poly = intersect(poly_roi,neighbours_polys{i}); 
            footp_area = area(neighbours_polys{i});
            neigh_coverage(i) = area(neigh_cov_poly)/footp_area;   
        catch
            neigh_coverage(i) = 0;
        end    
    end 

end