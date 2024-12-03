% 3D limb projection

% Input data
target = 'GANYMEDE';
obs    = 'JUICE';
et     = 1089532300.1838222;

% Parameters for cspice_limbpt function
flag = false; % assume the target area is visible from the instrument
method = 'TANGENT/ELLIPSOID';
[~, targetframe, ~] = cspice_cnmfrm(target); % body-fixed frame
abcorr = 'XLT+S';
corloc = 'CENTER';
refvec = [0; 0; 1]; % first of the sequence of cutting half-planes
ncuts  = 1e3; % number of cutting half-planes
delrol = cspice_twopi() / ncuts; % angular step by which to roll the
% cutting half-planes about the observer-target vector
schstp = 1.0d-4; % search angular step size
soltol = 1.0d-7; % solution convergence tolerance

% Limb calculation with cspice_limbpt function
[~, limb, ~, ~] = cspice_limbpt(method, target, et, targetframe, abcorr, ...
    corloc, obs, refvec, delrol, ncuts, schstp, soltol, ncuts); % limb
% points expressed in targetframe ref frame

% Variables definition
radii = cspice_bodvrd(target,'RADII',3);

% If this is the first time that the function is called within the
% simulation, then we have to initialize the figure
figure;
axes;
set(gcf, 'units', 'normalized', 'OuterPosition',...
    [0.0052,0.3139,0.4201,0.6532]);
hold on; grid minor; axis equal; box on; grid on;
set(gca,'Color','k','GridColor','w','MinorGridColor','w', ...
    'XColor','k','YColor','k','ZColor','k','FontSize',15)

% Ellipsoidal model of the target
[x, y, z] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 100);
% Print topography map
globe = surf(x,y,-z, 'FaceColor', .5*[1 1 1], 'EdgeColor', 'none');

try
    % Topography map of the target
    image = strcat(lower(target),'-map.jpg');
    % Load target image for texture map
    c_text = imread(image);
    if length(size(c_text)) > 2 % reads color, not needed...
        c_text = c_text(:, :, 1);
    end
    set(globe,'FaceColor','texturemap','CData',repmat(c_text,1,1,3), ...
        'FaceAlpha',1,'EdgeColor','none');
catch
end

% Body-fixed reference frame
quiver3(0, 0, 0, 1.5*radii(1), 0, 0, 'c', 'linewidth', 2)
quiver3(0, 0, 0, 0, 1.5*radii(2), 0, 'g', 'linewidth', 2)
quiver3(0, 0, 0, 0, 0, 1.5*radii(3), 'r', 'linewidth', 2)

% Axes definition
ax = gca;
ax.XAxis.TickLabelColor = 'k';
ax.YAxis.TickLabelColor = 'k';
ax.ZAxis.TickLabelColor = 'k';
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('3D limb projection')

% Spacecraft position
scpos = cspice_spkpos(obs, et, targetframe, 'NONE', target);
plot3(scpos(1), scpos(2), scpos(3), 'p', 'MarkerFaceColor',...
    [255 146 4]/255, 'MarkerEdgeColor', [255 146 4]/255)

% Limb projection onto the body surface
x = limb(1, :);
y = limb(2, :);
z = limb(3, :);
for i=1:length(limb)
    plot3(x(i), y(i), z(i), '.b', 'MarkerSize', 4)
end
view(scpos);

% % Limit axes if necessary
% cxlim = xlim();
% cylim = ylim();
% czlim = zlim();
% if max(abs(cxlim)) > 10*radii(1), xlim([-10*radii(1) 10*radii(1)]); end
% if max(abs(cylim)) > 10*radii(2), ylim([-10*radii(2) 10*radii(2)]); end
% if max(abs(czlim)) > 10*radii(3), zlim([-10*radii(3) 10*radii(3)]); end