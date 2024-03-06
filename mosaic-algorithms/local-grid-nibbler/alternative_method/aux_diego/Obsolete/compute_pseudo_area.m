% In order to work with complex shapes, a boundary box is created around
% the ROI being it aligned with the footprint at the starting epoch.

function  areapoints = compute_pseudo_area(areapoints_real,center_area,theta,cut_roi)

    % Rotate the ROI aligning it with the footprint
    
    for i = 1:size(areapoints_real,1)
        % Area vertex vector
        vector_real(i,:) = areapoints_real(i,:) - center_area;
        % Area vertex rotated
        vector_rotated(i,:) = [vector_real(i,1)*cos(theta) - vector_real(i,2)*sin(theta),  vector_real(i,1)*sin(theta) + vector_real(i,2)*cos(theta)];
        areapoints_rotated(i,:) = center_area + vector_rotated(i,:);
    end    
    
    % Compute the boundary limits of the rotated area
    
    box_xlim = [min(areapoints_rotated(:,1)),max(areapoints_rotated(:,1))];
    box_ylim = [min(areapoints_rotated(:,2)),max(areapoints_rotated(:,2))];
    
    % Define the boundary box of the rotated area
    
    areapoints_vertexes = [  box_xlim(1)  ,  box_ylim(2);
                             box_xlim(2)  ,  box_ylim(2);
                             box_xlim(2)  ,  box_ylim(1);
                             box_xlim(1)  ,  box_ylim(1);
                             box_xlim(1)  ,  box_ylim(2)];

    % Rotate the boudary box to the original orientation

    for i = 1:size(areapoints_vertexes,1)
        vector_vertex(i,:) = areapoints_vertexes(i,:) - center_area;
        % Area vertex rotated
        vector_vertex_rotated(i,:) = [vector_vertex(i,1)*cos(-theta) - vector_vertex(i,2)*sin(-theta),  vector_vertex(i,1)*sin(-theta) + vector_vertex(i,2)*cos(-theta)];
        areapoints(i,:) = center_area + vector_vertex_rotated(i,:);
    end


    if cut_roi == 1
        k = 1;
        pseudo_roi_old = areapoints;     
        for i = 1:size(pseudo_roi_old,1)
            if pseudo_roi_old(i,1)>180
                pseudo_roi_old(i,1) = pseudo_roi_old(i,1) - 360;
            end    
        end     
        areapoints = [];  
        areapoints = adapt_roi_cutted(pseudo_roi_old);   
    end
    
end