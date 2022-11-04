function ax = mapPlot(filename)
scrsz = get(groot,'ScreenSize');
xtick = -180:60:180;
ytick = -90:30:90;

W = scrsz(3); L = scrsz(4);
figure('Position',[100 100 W*0.8 L*0.7]);
ax = axes;
title('Topography Map','FontSize',20)
hold on; grid on; box on;
I = flip(imread(filename));
imagesc([-180, 180], [-90, 90], I)
colormap('gray');
set(gca, 'XTick', xtick, 'YTick', ytick, 'Layer', 'top', 'GridColor', 'w', 'GridAlpha', 1, 'GridLineStyle', ':');
set(gca, 'FontSize', 15)
xlim([min(xtick) max(xtick)]);
ylim([min(ytick) max(ytick)]);
xlabel('Longitude [°]')
ylabel('Latitude [°]')

end