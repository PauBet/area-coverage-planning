clc; clear; close all;

inputkernels;
addpath('C:\Users\roure\Desktop\science-optimizer\RESSlib\');
initSPICEv(fullK(METAKR));

startTime = cspice_str2et('2011 SEP 29 23:10:00.000 TDB');
endTime   = cspice_str2et('2011 SEP 29 23:50:00.000 TDB');
step = 60;
et = startTime:step:endTime;
instName = 'DAWN_FC2';
scName = 'DAWN';
targetName = 'VESTA';
[~, targetFrame, ~] = cspice_cnmfrm(targetName);


figure
hold on; grid minor; axis equal;
radii = cspice_bodvrd('VESTA','RADII',3);
[xe, ye, ze] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),100);
surf(xe, ye, ze, 'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.50 0.50 0.50])
for t=1:length(et)
    P  = cspice_spkpos(scName, et(t), targetFrame, 'NONE', targetName);
    [~, instFrame, boresight, bounds] = cspice_getfov(cspice_bodn2c(instName), 4);
    rotationMatrix = cspice_pxform(instFrame, targetFrame, et(t));
    for i=1:length(bounds)
        rot = rotationMatrix*bounds(:,i);
        bbounds(:,i) = normalize(rot, 'norm') + normalize(P,'norm');
    end
    pyramid_vertex = bbounds*norm(P);
    pyramid_vertex(:,end+1) = P;

    if t>1
        set(p,'visible','off')
    end
    plot3(P(1),P(2),P(3),'pb')
    face= [2 3 5; 1 2 5; 3 4 5; 4 1 5];
    p = patch('Faces',face,'Vertices',pyramid_vertex','Facecolor',[0.66 0.85 0.41], 'FaceAlpha', 0.5);
    drawnow
    pause
end