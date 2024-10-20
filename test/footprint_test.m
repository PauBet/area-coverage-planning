clc; clear; close all;

inputkernels;
addpath('C:\Users\roure\Desktop\science-optimizer\RESSlib\');
initSPICEv(fullK(METAKR));

startTime = cspice_str2et('2011 SEP 29 23:10:00.000 TDB');
endTime   = cspice_str2et('2011 SEP 29 23:30:00.000 TDB');
step = 60;
et = startTime:step:endTime;
instName = 'DAWN_FC2';
observer = 'DAWN';
targetName = 'VESTA';

mapPlot('vesta-map.png');
for i=1:1
    instFootprint(instName, targetName, observer, et(i));
end

function footprint = instFootprint(instName, targetName, scName, t)

[~, instFrame, boresight, bounds] = cspice_getfov(cspice_bodn2c(instName), 4);
method = 'ELLIPSOID';
[~, targetFrame, ~] = cspice_cnmfrm(targetName);
abcorr = 'NONE';
[dobs, ~] = cspice_spkpos(scName,t,targetFrame,abcorr,targetName);

[surfPoint, ~, ~, ~] = cspice_sincpt(method, targetName, t, targetFrame, abcorr, scName, instFrame, boresight);
[surfPoint1, ~, ~, ~] = cspice_sincpt(method, targetName, t, targetFrame, abcorr, scName, instFrame, bounds(:,1));
[surfPoint2, ~, ~, ~] = cspice_sincpt(method, targetName, t, targetFrame, abcorr, scName, instFrame, bounds(:,2));
[surfPoint3, ~, ~, ~] = cspice_sincpt(method, targetName, t, targetFrame, abcorr, scName, instFrame, bounds(:,3));
[surfPoint4, ~, ~, ~] = cspice_sincpt(method, targetName, t, targetFrame, abcorr, scName, instFrame, bounds(:,4));

footprint = [surfPoint1, surfPoint2, surfPoint3, surfPoint4];

fp.ox = surfPoint(1);
fp.oy = surfPoint(2);
fp.oz = surfPoint(3);

mid = zeros(3,length(bounds));
for i=1:length(bounds)
    if i~=length(bounds)
        mid(:,i) = .5*(bounds(:,i) + bounds(:,i+1));
    else
        mid(:,i) = .5*(bounds(:,i) + bounds(:,1));
    end
end

fp.sizex = .5*norm((mid(:,2) - mid(:,4)));
fp.sizey = .5*norm((mid(:,1) - mid(:,3)));
fp.overlapx = 0.1;
fp.overlapy = 0.8;

[~, lon, lat] = cspice_reclat(surfPoint);
plot(lon, lat,'y^')
hold on;
[dobs, ~] = cspice_subpnt('INTERCEPT/ELLIPSOID',targetName,t,targetFrame,abcorr,scName);
[~, lon, lat] = cspice_reclat(dobs);
plot(lon, lat, 'm^')

[~, lon1, lat1] = cspice_reclat(surfPoint1);
[~, lon2, lat2] = cspice_reclat(surfPoint2);
[~, lon3, lat3] = cspice_reclat(surfPoint3);
[~, lon4, lat4] = cspice_reclat(surfPoint4);
areapoints = [ lon1 lat1;
               lon2 lat2;
               lon3 lat3;
               lon4 lat4];
areapoints = areapoints.*cspice_dpr;
targetArea = polyshape(areapoints);
plot(targetArea, 'FaceColor', [0.93,0.69,0.13])
legend('Target area')

end