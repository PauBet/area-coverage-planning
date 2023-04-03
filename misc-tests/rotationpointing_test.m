clc; clear; close all;

inputkernels;
addpath('C:\Users\roure\Desktop\science-optimizer\RESSlib\');
initSPICEv(fullK(METAKR));

startTime = cspice_str2et('2011 SEP 29 23:10:00.000 TDB');
endTime   = cspice_str2et('2011 SEP 29 23:30:00.000 TDB');
step = 60;
et = startTime:step:endTime;
et = startTime;
instName = 'DAWN_FC2';
scName = 'DAWN';
targetName = 'VESTA';
targetFrame = 'IAU_VESTA';
[~, instFrame, boresight, bounds] = cspice_getfov(cspice_bodn2c(instName), 4);

% initial spacecraft clocktime
sclock = cspice_sce2c(cspice_bodn2c(scName),et(1));

%function to obtain the matrix with the initial orientation of the
%instrument
%cmatfunc = @() cspice_ckgpav(instName, sclock, 5.*256.,target_fixed);
cmatfunc = cspice_ckgpav(cspice_bodn2c(instName), sclock, 5.*256.,instFrame);

%% target point
px = -130*cspice_rpd;
py = 30*cspice_rpd;
latPoint = [mean(cspice_bodvrd('VESTA','RADII',3)), px, py];
recPoint = cspice_latrec(latPoint(1), latPoint(2), latPoint(3));
instPos  = cspice_spkpos(scName, et, targetFrame, 'NONE', targetName);
v2 = recPoint - instPos;
rotationMatrix = cspice_pxform(targetFrame, instFrame, et);
v2 =  rotationMatrix*v2;
rotAxis = normalize(cross(boresight, v2), 'norm');
angle = cspice_vsep(boresight, v2);
pointingRotation = cspice_axisar(rotAxis, angle);
newPointing = pointingRotation*cmatfunc;

for i=1:length(bounds)
    newbounds(:,i) = pointingRotation*bounds(:,i);
end

method = 'ELLIPSOID';
abcorr = 'NONE';
[oldsurfPoints(:,1), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, bounds(:,1));
[oldsurfPoints(:,2), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, bounds(:,2));
[oldsurfPoints(:,3), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, bounds(:,3));
[oldsurfPoints(:,4), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, bounds(:,4));
[oldsurfPoints(:,5), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, boresight);
for i=1:length(oldsurfPoints)
    [~, olatSurfPoints(i,1), olatSurfPoints(i,2)] = cspice_reclat(oldsurfPoints(:,i));
    olatSurfPoints(i,:) = olatSurfPoints(i,:)*cspice_dpr;
end

[surfPoints(:,1), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, newbounds(:,1));
[surfPoints(:,2), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, newbounds(:,2));
[surfPoints(:,3), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, newbounds(:,3));
[surfPoints(:,4), ~, ~, ~] = cspice_sincpt(method, targetName, et, targetFrame, abcorr, scName, instFrame, newbounds(:,4));

for i=1:length(surfPoints)
    [~, latSurfPoints(i,1), latSurfPoints(i,2)] = cspice_reclat(surfPoints(:,i));
    latSurfPoints(i,:) = latSurfPoints(i,:)*cspice_dpr;
end

%%
mapPlot('vesta-map.png')
plot(px*cspice_dpr, py*cspice_dpr, 'y^', 'MarkerSize', 8);
hold on; box on;
for i=1:length(latSurfPoints)
    plot(latSurfPoints(i,1), latSurfPoints(i,2),'m^','MarkerSize', 8)
end
for i=1:length(oldsurfPoints)
    plot(olatSurfPoints(i,1), olatSurfPoints(i,2),'c^','MarkerSize', 8)
end
