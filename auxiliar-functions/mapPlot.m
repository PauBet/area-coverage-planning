function ax = mapPlot(filename)
scrsz = get(groot,'ScreenSize');
xtick = -180:60:180;
ytick = -90:30:90;

W = scrsz(3); L = scrsz(4);
figure('Position',[0.4307*W 0.3144*L W*0.6 L*0.5]);
ax = axes;
title('Coverage Map')
hold on; grid on; box on; axis equal;
I = flip(imread(filename));
imagesc([-180, 180], [-90, 90], I)
colormap('gray');
set(gca, 'XTick', xtick, 'YTick', ytick, 'Layer', 'top', ...
    'GridColor', 'w', 'GridAlpha', 1, 'GridLineStyle', ':');
set(gca, 'FontSize', 20)
xlim([min(xtick) max(xtick)]);
ylim([min(ytick) max(ytick)]);
xlabel('Planetocentric longitude [°]')
ylabel('Planetocentric latitude [°]')

end