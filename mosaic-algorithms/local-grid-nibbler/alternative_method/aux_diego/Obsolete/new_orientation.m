%% NEW POINTING MATRIX

function [rotation] = new_orientation(et,zerorotation,targetcoord)

Europarad = cspice_bodvrd('EUROPA','RADII',3);

% initial target: original target

zinitial = zerorotation.'*[0; 0; 1]; 
xinitial = zerorotation.'*[1; 0; 0];
yinitial = zerorotation.'*[0; 1; 0];

% Conversion of the first target from lonlat to rectangular coordinates

midpointrec = cspice_latrec(mean(Europarad),targetcoord(1)*(pi/180),targetcoord(2)*(pi/180)); %Center target in IAU_EUROPA rectangular
 
% Obtetion of the vector between the first target and the spacecraft
% position in the IAU_Europa reference frame

[dgali , ~] = cspice_spkpos('-77',et,'IAU_EUROPA','NONE','EUROPA');

zvectorIAU = midpointrec-dgali;

%Distance to the target 

distance = norm(zvectorIAU);

%Normalizing the previously obtained vector

zvectorIAU = zvectorIAU/(norm(zvectorIAU));

% Rotation axis

rotaxis_0 = cross(zinitial,zvectorIAU);
rotaxis_0 = rotaxis_0/(norm(rotaxis_0));

% angle

phi = acos(dot(zinitial,zvectorIAU));

% Rotation matrix (Rodrigues rotation formula)

w_0 = [0            -rotaxis_0(3)  rotaxis_0(2);
       rotaxis_0(3)  0            -rotaxis_0(1);
      -rotaxis_0(2)  rotaxis_0(1)  0];

rotation_0 = eye(3) + sin(phi)*w_0 + (2*sin(phi/2)^2)*(w_0)^2;

%final x and y instrument axis

xvectorIAU = rotation_0*xinitial;
yvectorIAU = rotation_0*yinitial;

% tranformation matrix new instrument frame

rotation = [xvectorIAU(1) yvectorIAU(1) zvectorIAU(1); xvectorIAU(2) yvectorIAU(2) zvectorIAU(2); xvectorIAU(3) yvectorIAU(3) zvectorIAU(3)];


end