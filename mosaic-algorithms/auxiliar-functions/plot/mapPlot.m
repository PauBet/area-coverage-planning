function ax = mapPlot(filename)

% Future check: filename does not exist!

scrsz = get(groot,'ScreenSize');
xtick = -180:45:180;
ytick = -90:30:90;
% x tick label
xtickstr = string();
for i=1:length(xtick)
    if xtick(i) < 0 && xtick(i) > -180
        xtickstr(i) = strcat(num2str(-xtick(i)), 'º', 'W');
    elseif xtick(i) > 0 && xtick(i) < 180
        xtickstr(i) = strcat(num2str(xtick(i)), 'º', 'E');
    else
        xtickstr(i) = strcat(num2str(abs(xtick(i))), 'º');
    end
end

% y tick label
ytickstr = string();
for i=1:length(ytick)
    if ytick(i) < 0
        ytickstr(i) = strcat(num2str(-ytick(i)), 'º', 'S');
    elseif ytick(i) > 0
        ytickstr(i) = strcat(num2str(ytick(i)), 'º', 'N');
    else
        ytickstr(i) = strcat(num2str(ytick(i)), 'º');
    end
end

% figure
W = scrsz(3); L = scrsz(4);
figure('Position',[0.4307*W 0.3144*L W*0.6 L*0.5]);
ax = axes;
hold on; grid on; box on; axis equal;
map = flip(imread(filename));
%im = image(map);
%I = flip(imread(filename));
im = imagesc([-180, 180], [-90, 90], map);
im.AlphaData = 0.8;
colormap('gray');
set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
    'YTickLabel', ytickstr, 'Layer', 'top', 'GridColor', 'w', ...
    'GridAlpha', 1, 'GridLineStyle', ':');
set(gca, 'FontSize', 20)
xlim([min(xtick) max(xtick)]);
ylim([min(ytick) max(ytick)]);
xlabel('Planetocentric longitude')
ylabel('Planetocentric latitude')

end