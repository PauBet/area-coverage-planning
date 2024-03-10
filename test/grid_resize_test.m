close all;
% grid resize test

N = 11;
x = linspace(1, 2, N);
y = linspace(1, 2, N);
map = cell(N, N);
for i=1:N
    for j=1:N
        map{i, j} = [x(j), y(i)];
    end
end

figure
hold on; grid minor;
for i=1:N
    for j=1:N
        plot(map{i, j}(1), map{i, j}(2), 'k*')
    end
end

step = 0.1;
delta = 0.1;

for i=1:N
    for j=1:N
        if i == 1
            if j ~= 1
                map{i, j}(1) =  map{i, j-1}(1) + step + delta;
            end
        elseif j == 1
            if i ~= 1
                map{i, j}(2) = map{i-1, j}(2) + step + delta;
            end
        else
            map{i, j}(1) = map{i, j-1}(1) + step + delta;
            map{i, j}(2) = map{i-1, j}(2) + step + delta;
        end
        
    end
end

hold on; grid minor;
for i=1:N
    for j=1:N
        plot(map{i, j}(1), map{i, j}(2), 'r*')
    end
end

% Orient map (30º)
map = cell(N, N);
for i=1:N
    for j=1:N
        map{i, j} = [x(j), y(i)];
    end
end

angle = deg2rad(30);
rotmat = [cos(angle)   -sin(angle);
          sin(angle)   cos(angle)];
% matrixGrid directions x and y
dirx = rotmat(1, :);
diry = rotmat(2, :);
cx = 1.5; cy = 1.5;
orientedMap  = cell(N, N);
for i=1:N
    for j=1:N
        orientedMap{i, j} = [cx, cy]' + rotmat*(map{i, j}' - ...
            [cx, cy]');
    end
end

figure
hold on; grid minor;
for i=1:N
    for j=1:N
        plot(orientedMap{i, j}(1), orientedMap{i, j}(2), 'k*')
    end
end

for i=1:N
    for j=1:N
        if i == 1
            if j ~= 1
                orientedMap{i, j}(1) =  orientedMap{i, j-1}(1) + (step + delta)*dirx(1);
                orientedMap{i, j}(2) =  orientedMap{i, j-1}(2) + (step + delta)*diry(1);
            end
        elseif j == 1
            if i ~= 1
                orientedMap{i, j}(1) = orientedMap{i-1, j}(1) + (step+delta)*dirx(2);
                orientedMap{i, j}(2) = orientedMap{i-1, j}(2) + (step+delta)*diry(2);
            end
        else
            aux1(1) = orientedMap{i, j-1}(1) + (step + delta)*dirx(1);
            aux1(2) = orientedMap{i, j-1}(2) + (step + delta)*diry(1);
            % aux1 == aux2;
            % aux2(1) = orientedMap{i-1, j}(1) + (step+delta)*dirx(2);
            % aux2(2) = orientedMap{i-1, j}(2) + (step+delta)*diry(2);
            orientedMap{i, j} = aux1;
        end
        
    end
end

hold on; grid minor;
for i=1:N
    for j=1:N
        plot(orientedMap{i, j}(1), orientedMap{i, j}(2), 'r*')
    end
end

% Irregular grid
% figure
% hold on; grid minor; axis equal;
% set(gca, 'xlim', [0 4], 'ylim', [0 4], 'fontsize', 20)
% ax = gca;
% roipoly = drawpolygon(ax);
% targetArea = roipoly.Position;

clear orientedArea aux;

targetArea = [1.1095    2.4905;
    2.0088    3.2730;
    3.1066    2.2569;
    3.0015    1.4511;
    2.0905    0.5518;
    1.1562    1.3109;
    1.7752    1.9299];

[cx, cy] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
gamma = [cx, cy];

stepx = 0.5; stepy = 0.3; deltax = 0.01; deltay = 0.02;
fpref.sizex = stepx; fpref.sizey = stepy; fpref.angle = 30;
ovlapx = 10; ovlapy = 10;
[matrixGrid, dirx, diry] = grid2D(fpref, ovlapx, ovlapy, ...
    gamma, targetArea);

figure
plot(polyshape(targetArea(:, 1), targetArea(:, 2)))
hold on; box on;
for i=1:size(matrixGrid, 1)
    for j=1:size(matrixGrid, 2)
        if ~isempty(matrixGrid{i, j})
            plot(matrixGrid{i, j}(1), matrixGrid{i, j}(2), 'k*')
        end
    end
end

orientedMap = matrixGrid;
Nx = size(orientedMap, 1);
Ny = size(orientedMap, 2);
refel = orientedMap{1, 2};

for i=1:Nx
    for j=1:Ny
        if ~isempty(orientedMap{i, j})
            if norm(orientedMap{i, j} - refel) < 1e-5
                indrow = i; indcol = j;
                flag = true;
                break;
            end
        end
    end
    if flag, break; end
end

for i=1:Nx
    for j=1:Ny
        deltarow = i - indrow;
        deltacol = j - indcol;
        if ~isempty(orientedMap{i ,j})
            orientedMap{i, j}(1) =  orientedMap{indrow, indcol}(1) + deltacol*(stepx + deltax)*dirx(1) - deltarow*(stepy + deltay)*diry(1);
            orientedMap{i, j}(2) =  orientedMap{indrow, indcol}(2) + deltacol*(stepx + deltax)*dirx(2) - deltarow*(stepy + deltay)*diry(2);
            plot(orientedMap{i, j}(1), orientedMap{i, j}(2), 'r*')
        end
    end
end

% particular case angle = 180º
fpref.angle = 180;
olapx = ovlapx*fpref.sizex/100; olapy = ovlapy*fpref.sizey/100;
[matrixGrid, dirx, diry] = grid2D(fpref, ovlapx, ovlapy, ...
    gamma, targetArea);

figure
hold on;
plot(polyshape(targetArea(:, 1), targetArea(:, 2)));
for i=1:Nx
    for j=1:Ny
        if ~isempty(matrixGrid{i, j})
            plot(matrixGrid{i, j}(1), matrixGrid{i, j}(2), 'bo')
        end
    end
end

map = grid2Dresize(gamma', matrixGrid, {}, fpref.sizex - olapx, fpref.sizey - olapy, dirx, diry, ...
    fpref.sizex - olapx, fpref.sizey - olapy);

hold on;
for i=1:Nx
    for j=1:Ny
        if ~isempty(matrixGrid{i, j})
            plot(map{i, j}(1), map{i, j}(2), 'r*')
        end
    end
end

