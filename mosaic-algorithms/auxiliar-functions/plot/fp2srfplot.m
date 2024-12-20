function fp2srfplot(fp, videosave)
% This figure shows the 3D FOV projection of the instrument over the target's
% surface, modeled as a triaxial ellipsoid
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% Version:      2
% Last update:  12/2023
% 
% Usage:        fp2srfplot(fp, videosave)
%
% Inputs:
%   > fp:       struct containing the parameters that define the
%               footprint. See function footprint for further information

% Variables definition
persistent counter; % function call counter
persistent ax; % figure axes definition as a persistent variable to plot
% more than one footprint
persistent video;
[~, targetframe, ~] = cspice_cnmfrm(fp.target);
radii = cspice_bodvrd(fp.target,'RADII',3);

% This figure may be additive, i.e., could be iteratively called in order
% to print different footprints
if isempty(counter) 
    counter = 0;
    if videosave
        video = VideoWriter('footprint_3dprojection');
        video.FrameRate = 1;
        open(video);
    end
else 
    counter = counter + 1;
    if videosave, open(video); end
end

if ~counter
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
        image = strcat(lower(fp.target),'-map.png');
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
    title('Footprint projection')
else
    hold on;
end

% Instrument (spacecraft) position
instpos = cspice_spkpos(fp.sc, fp.t, targetframe, 'NONE', fp.target);
plot3(instpos(1), instpos(2), instpos(3), 'p', 'MarkerFaceColor',...
    [255 146 4]/255, 'MarkerEdgeColor', [255 146 4]/255)

% Perspective projection of the FOV onto the body surface (pyramid)
fp.fovbounds = normalize(fp.fovbounds, 'norm') + normalize(instpos, 'norm');
pyramid_vertex = norm(instpos)*fp.fovbounds;
pyramid_vertex(:,end+1) = instpos;
face = [2 3 5; 1 2 5; 3 4 5; 4 1 5];
p = patch('Faces',face,'Vertices',pyramid_vertex','Facecolor',...
     [0.66 0.85 0.41], 'FaceAlpha', 0.5);
hold on;

% Footprint projection onto the body surface
x = fp.recVertices(:, 1);
y = fp.recVertices(:, 2);
z = fp.recVertices(:, 3);
for i=1:length(fp.recVertices)
    if ~isnan(fp.recVertices)
        plot3(x(i), y(i), z(i), '.b', 'MarkerSize', 4)
    end
end
view(instpos);

% Limit axes if necessary
cxlim = xlim();
cylim = ylim();
czlim = zlim();
if max(abs(cxlim)) > 10*radii(1), xlim([-10*radii(1) 10*radii(1)]); end
if max(abs(cylim)) > 10*radii(2), ylim([-10*radii(2) 10*radii(2)]); end
if max(abs(czlim)) > 10*radii(3), zlim([-10*radii(3) 10*radii(3)]); end

% Save in video
if videosave
    writeVideo(video, getframe(gcf));
    set(p, 'Visible', 'off')
end


% Close video
close(video);

end