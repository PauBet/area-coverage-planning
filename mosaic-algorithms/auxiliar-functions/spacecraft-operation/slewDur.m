function [tdur, slewOk] = slewDur(p1, p2, t, inst, target, sc, slew_rate)

[~, ~, R1, ~, ~] = instpointing(inst, target, sc, t, p1(1), p1(2));
[~, ~, R2, ~, ~] = instpointing(inst, target, sc, t, p2(1), p2(2));

Rdelta = transpose(R1)*R2;
angle = acos((trace(Rdelta) - 1)/2);
tdur = angle/slew_rate;

end