% In some cases the footprint will be rotated, considering that the
% rotation cannot be deliberately modified, the mosaic generation will be
% addapted to this rotation. For this reason it is necessary to compute the
% principal axes of the footprint and its rotation with respect to the
% standard 2d reference frame

function [theta, x_foot, y_foot] = foot_axes(zerofootprint)

    % Compute the mean points on the 4 sides of the footprint
    
    mean_p1 = 0.5*[zerofootprint(1,54)+zerofootprint(1,154), zerofootprint(2,54)+zerofootprint(2,154)];
    mean_p2 = 0.5*[zerofootprint(1,154)+zerofootprint(1,254), zerofootprint(2,154)+zerofootprint(2,254)];
    mean_p3 = 0.5*[zerofootprint(1,254)+zerofootprint(1,354), zerofootprint(2,254)+zerofootprint(2,354)];
    mean_p4 = 0.5*[zerofootprint(1,354)+zerofootprint(1,54), zerofootprint(2,354)+zerofootprint(2,54)];
    
    % Estimate the footprint principal axes
    
    x1_foot = (mean_p1-mean_p3)/norm(mean_p1-mean_p3);
    x2_foot = (mean_p2-mean_p4)/norm(mean_p2-mean_p4);
    
    % The outputs are defined following these criteria in order to ensure
    % efectiveness
    
    if abs(x2_foot(1))>abs(x1_foot(1))
        if x2_foot(1)<0
            x_foot = -1*x2_foot;
            y_foot = x1_foot;
        else
            x_foot = x2_foot;
            y_foot = x1_foot;
        end    
    else
        if x1_foot(1)<0
            x_foot = -1*x1_foot;
            y_foot = x2_foot;
        else
            x_foot = x1_foot;
            y_foot = x2_foot;
        end  
    end
    
    theta = 1.1*acos(dot([1,0],x_foot));

end    