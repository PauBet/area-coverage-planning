a = [2, 3];
b = [-3, 5];

figure
plot([a(1), b(1)], [a(2), b(2)], 'linewidth', 1)
hold on; box on; grid minor;
plot(a(1), a(2), 'k*')
plot(b(1), b(2), 'k*')
set(gca, 'fontsize', 15)

alpha = alphaslope(a, b);
delta = [cos(alpha), sin(alpha)];

plot(a(1) + delta(1), a(2) + delta(2), 'r*')

function alpha = alphaslope(p1, p2)
% This function outputs the slope of the a line defined by two points
alpha = atan2(p2(2) - p1(2), p2(1) - p1(1));
end