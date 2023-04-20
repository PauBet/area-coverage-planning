function plot_area_and_footprints(startTime,real_targets,roi,inst,target_fixed,footprint_func)

    zero_orientation = cspice_pxform(inst,target_fixed,startTime);

    %Plot the final footprints

    close all
    fill(roi(:,1),roi(:,2),'r','EdgeColor','#808080','FaceColor','#808080')
    hold on
    for i = 1:size(real_targets,1)
        if i == 1
            target_orientation = new_orientation(real_targets(i,3),zero_orientation.',real_targets(i,1:2));
        else
            target_orientation = new_orientation(real_targets(i,3),target_orientation.',real_targets(i,1:2));
        end    
        target_footprint = footprint_func(real_targets(i,3),target_orientation,200);
        plot(target_footprint(1,:),target_footprint(2,:),'LineWidth',1.5,'Color','#FF0000')
    end
end