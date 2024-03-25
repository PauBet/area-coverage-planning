function fp_coverage = get_fp_coverage(poly_roi,fp_poly)
% Computes the ratio of coverage of a footprint with respect to the ROI, it uses 
% polyshape objects to compute it. 
%
% Programmers: Diego Andía (UPC/ESEIAAT)
% Date:        08/2023
% Version:      1
% Last update: 03/2024
%
% Usage:       fp_coverage = get_fp_coverage(poly_roi,fp_poly)
%
% Inputs: 
%    > poly_roi: polyshape object containing the (remaining) ROI
%    > fp_poly: polyshape object conatining the footprint
% Outputs:
%    > fp_coverage: Footprint's ratio of ROI coverage 
    
    try
       intersec = intersect(poly_roi,fp_poly); 
       intersec_area = area(intersec);
       fp_coverage = intersec_area/area(fp_poly);   
    catch
       fp_coverage = 0;
    end    

end