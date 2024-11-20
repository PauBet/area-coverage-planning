function fp_corners = get_footprint_corners(footprint)
    
    n = 1;

    for i = 2:length(footprint)
        if norm(footprint(:,i)-footprint(:,i-1)) == 0
            fp_corners(:,n) = footprint(:,i);
            n = n + 1;
        end    
    end    

end
