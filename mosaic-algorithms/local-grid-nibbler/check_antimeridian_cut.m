function [roi, cut_roi] = check_antimeridian_cut(roi)
% Given a ROI obtains its intersection points with the antimeridian line
% and splits the ROI into two ROIs contained in the range of longitudes -180 +180
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date:        08/2023
% Version:      1
% Last update: 03/2024
%
% Usage:       [roi, roi_real, cut_roi] = check_antimeridian_cut(roi)
%
%
% Inputs: 
%    > roi: Coordinates in lon-lat of the region of interest
% Outputs:
%    > roi: Coordinates in lon-lat of the region of interest after
%           applying the required adaptation                                                    
%    > cut_roi: Flag that indicates if the roi has been splitted due 
%           to antimeridian cut (1 -> True)

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
        roi = adapt_roi_cutted(roi_old);
    end

end