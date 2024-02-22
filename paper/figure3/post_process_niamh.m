% Zoom in
xlim([120 180])
ylim([0 40])
xtick = 120:15:180;
ytick = 0:10:40;

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
set(gca, 'XTick', xtick, 'YTick', ytick, 'XTickLabel', xtickstr, ...
    'YTickLabel', ytickstr, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Layer', 'top', 'GridColor', 'w', 'GridAlpha', 1, 'GridLineStyle', ':');

% Save figure [PDF]
figpath = '.';
set(gcf, 'Units', 'inches', 'Position', [3,3,7.3,4.9]);
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);
roiname = char(lower(roistruct(1).name));
roiname(isspace(roiname)) = [];
name = ['post_process_',roiname];
filename = fullfile(figpath, [roiname,'_',mosaic]);
print(gcf, filename, '-dpdf', '-r1200')