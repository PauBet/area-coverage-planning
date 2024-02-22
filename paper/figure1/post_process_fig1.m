% Zoom in
xlim([-30 60])
ylim([-15 35])
xtick = -30:15:60;
ytick = -15:10:35;

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