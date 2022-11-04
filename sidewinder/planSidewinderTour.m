function tour = planSidewinderTour(closestSide, vertices, fprint0, ...
    gamma, olapx, olapy)
% This function discretizes and computes the coverage path of a certain
% ROI in order to build a mosaic image, adapted from [1].
%
% Programmers:  Paula Betriu (UPC/ESEIAAT)
% Date:         09/2022
% 
% Usage:        tour = planSidewinderTour(closestSide, vertices,...
%               fprint0, gamma)
%
% Inputs:
%   > closestSide:  string defining the closest side of the polygon (see
%                   vertices) to the spacecraft ground track position
%   > vertices:     matrix containing the vertices of the ROI polygon. The
%                   vertex points are expressed in 2D. 
%       # vertices(:,1) correspond to the x values of the vertices
%       # vertices(:,2) correspond to the y values of the vertices
%   > fprint0:      struct containing the footprint parameters that are
%                   going to be used to define the grid discretization
%   > gamma:        origin of the grid, in latitudinal coordinates, in deg
%   > olapx:        grid footprint overlap in the x direction (longitude),
%                   in deg
%   > olapy:        grid footprint overlap in the y direction (latitude),
%                   in deg
% 
% Outputs:
%   > tour:         cell matrix of the successive planned observations.
%                   Each observation is defined by the instrument boresight
%                   projection onto the body surface, in latitudinal
%                   coordinates [lon lat], in deg
%
% [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G., 
% Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis 
% steerable, 2D framing sensors.

% Pre-allocate variables
tour = {}; % list of planned observations
bearing = true; % boolean that will alternate coverage path directions

% 2D grid discretization
%grid = grid2D(fprint0.sizex, fprint0.sizey, olapx, olapy, gamma, vertices);
[grid, vlon, vlat] = grid2D(fprint0.sizex, fprint0.sizey, olapx, ...
    olapy, gamma, vertices);

% The origin of the coverage path depends on the spacecraft ground track
% position
switch closestSide
    case {'north','south'} % Horizontal sweep
        for i=1:size(grid,1)
            % Sweep across latitude
            y = vlat(length(vlat) - i + 1);
            %y = grid{i,1}(2); % Highest to lowest latitude
            if isequal(closestSide, 'north')
                %y = grid{size(grid,1) + 1 - i, 1}(2); % Lowest to highest
                % latitude
                y = vlat(i);
            end
            for j=1:size(grid,2)
                % Sweep across longitude
                x = grid{i,size(grid,2) + 1 - j}(1); % Highest to lowest
                % longitude
                if bearing
                    x = grid{i,j}(1); % Lowest to highest longitude
                end
                % The boundary box may not correspond to the target area,
                % therefore the grid matrix -sized as the boundary box- may
                % contain some points that are outside of the target
                % polygon, saved as [NaN NaN]. See grid2D for further
                % information
                if ~isnan(x)
                    tour{end + 1} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each 
            % row sweeping, i.e. left (highest lon) to right (lowest lon) 
            % or vice versa
        end
    case {'east','west'} % Vertical sweep
        for i=1:size(grid,2)
            % Sweep across longitude
            x = vlon(i);
            %x = grid{1,size(grid,2) + 1 - i}(1); % Lowest to highest
            % latitude
            if isequal(closestSide, 'east')
                %x = grid{1,i}(1); % Highest to lowest latitude
                x = vlon(length(vlon) + i - 1);
            end
            % The boundary box may not correspond to the target area,
            % therefore the grid matrix -sized as the boundary box- may
            % contain some points that are outside of the target
            % polygon, saved as [NaN NaN]. See grid2D for further
            % information
            for j=1:size(grid,1)
                % Sweep across latitude
                y = grid{j,i}(2); % Highest to lowest latitude
                if bearing
                    y = grid{size(grid,1) + 1 - j, i}(2); % Lowest to
                    % highest latitude
                end
                if ~isnan(y)
                    tour{end + 1} = [x y]; % Save it in the coverage tour
                end
            end
            bearing = not(bearing); % Switch coverage direction after each
            % column sweeping, i.e. north (highest lat) to south (lowest
            % lat) or vice versa
        end
end

%% Obsolete
% switch closestSide
%     case {'north','south'}
%         for i=1:length(yi)
%             y = ymax - yi(i) + ymin;
%             if isequal(closestSide, 'north')
%                 y = yi(i);
%             end
%             for j=1:length(xi)
%                 x = xmax - xi(j) + xmin;
%                 if bearing
%                     x = xi(j);
%                 end
%                 in = inPolygonCheck(x, y, vertices(:,1), vertices(:,2), [bbox.minlon bbox.maxlon], [bbox.minlat bbox.maxlat]);
%                 if in
%                     tour{end + 1} = [x y];
%                 end
%             end
%             bearing = not(bearing);
%         end
%     case {'east','west'}
%         for i=1:length(xi)
%             x = xmax - xi(i) + xmin;
%             if isequal(closestSide, 'east')
%                 x = xi(i);
%             end
%             for j=1:length(yi)
%                 y = ymax - yi(j) + ymin;
%                 if bearing
%                     y = yi(j);
%                 end
%                 in = inPolygonCheck(x, y, vertices(:,1), vertices(:,2), [bbox.minlon bbox.maxlon], [bbox.minlat bbox.maxlat]);
%                 if in
%                     tour{end + 1} = [x y];
%                 end
%             end
%             bearing = not(bearing);
%         end
% end

end