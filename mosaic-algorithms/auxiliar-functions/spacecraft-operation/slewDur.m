function [tdur, slewOk] = slewDur(v1, v2, slew_rate)

angle = cspice_vsep(v1, v2);
tdur = angle/slew_rate;

end