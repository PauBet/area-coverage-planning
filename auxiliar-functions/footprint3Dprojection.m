function footprint3Dprojection(fp, videosave)
% This figure shows the FOV projection of the instrument over the target's
% surface, modeled as a triaxial ellipsoid
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         10/2022
% 
% Usage:        footprint3Dprojection(fp)
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
    set(gcf,'units','normalized','OuterPosition',...
        [0.0052,0.3139,0.4201,0.6532]);
    hold on; grid minor; axis equal; box on; grid on;
    set(gca,'Color','k','GridColor','w','MinorGridColor','w', ...
        'XColor','k','YColor','k','ZColor','k','FontSize',15)

    % Ellipsoidal model of the target
    radii = cspice_bodvrd(fp.target,'RADII',3);
    [x, y, z] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 100);

    % Topography map of the target
    image = strcat(lower(fp.target),'-map.png');
    % Load target image for texture map
    c_text = imread(image);
    if length(size(c_text)) > 2 % reads color, not needed...
        c_text = c_text(:, :, 1);
    end
    % Print topography map
    globe = surf(x,y,-z, 'FaceColor', .5*[1 1 1], 'EdgeColor', .5*[1 1 1]);
    set(globe,'FaceColor','texturemap','CData',repmat(c_text,1,1,3), ...
        'FaceAlpha',1,'EdgeColor','none');

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
plot3(ax, instpos(1),instpos(2),instpos(3),'p','MarkerFaceColor',...
    [255 146 4]/255, 'MarkerEdgeColor', [255 146 4]/255)

% Perspective projection of the FOV onto the body surface (pyramid)
fp.fovbounds = normalize(fp.fovbounds, 'norm') + normalize(instpos, 'norm');
pyramid_vertex = norm(instpos)*fp.fovbounds;
pyramid_vertex(:,end+1) = instpos;
face = [2 3 5; 1 2 5; 3 4 5; 4 1 5];
p = patch(ax, 'Faces',face,'Vertices',pyramid_vertex','Facecolor',...
     [0.66 0.85 0.41], 'FaceAlpha', 0.5);
hold on;

% Footprint projection onto the body surface
for i=1:length(fp.bvertices)
    if ~isnan(fp.bvertices(i, :))
        surfPoint = cspice_srfrec(cspice_bodn2c(fp.target), ...
            fp.bvertices(i, 1)*cspice_rpd, fp.bvertices(i, 2)*cspice_rpd);
        plot3(ax, surfPoint(1), surfPoint(2), surfPoint(3),...
            '.b','MarkerSize',4)
    end
end
view(instpos);

% Save in video
if videosave, writeVideo(video, getframe(gcf)); end
set(p, 'Visible', 'off')

% Close video
close(video);

end