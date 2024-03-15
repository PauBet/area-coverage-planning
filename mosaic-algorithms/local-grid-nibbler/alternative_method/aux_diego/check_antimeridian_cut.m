function [roi, roi_real, cut_roi] = check_antimeridian_cut(roi)

    max_lon = max(roi(:,1)); min_lon = min(roi(:,1)); 

    cut_roi = 0;

    if abs(max_lon-min_lon)>200 
        cut_roi = 1;
    end    

    if cut_roi == 1        
        roi_old = roi;
        for i = 1:size(roi,1)
            if roi(i,1)<0
                roi(i,1) = roi(i,1) + 360;
            end    
        end         
        roi_real = roi;
        roi = [];
        roi = adapt_roi_cutted(roi_old);
    else
        roi_real = roi;
    end   

end