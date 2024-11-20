function [areapoints] = add_area_manual(n)

for i = 1:n
    areapoints(i,:) = ginput(1);
    plot(areapoints(i,1) , areapoints(i,2),'.','MarkerSize',8,'Color','blue')
    hold on
    axis manual
end

areapoints(end+1,:) = areapoints(1,:) ;
plot(areapoints(:,1) , areapoints(:,2),'LineWidth',2.5,'Color','blue')

end