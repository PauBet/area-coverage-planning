function grid_topo = inst2topo(grid, lon, lat, target, sc, inst, et)

% Pre-allocate
[~, ~, rotmat] = instpointing(inst, target, sc, et, lon, lat);
grid_topo = cell(size(grid));
method = 'ELLIPSOID';
[~, targetframe, ~] = cspice_cnmfrm(target); % target frame ID in SPICE
instpos  = cspice_spkpos(sc, et, targetframe, 'NONE', target); % rectangular
            % coordinates of the instrument in the body-fixed reference frame

% Convert grid into topographical coordinates
for i=1:size(grid, 1)
    for j=1:size(grid, 2)
        instp = grid{i, j};
        if ~isempty(instp) && ~any(isnan(instp), 'all')
            p = zeros(3, 1);
            p(1:2) = instp; p(3) = 1;
            p_body = rotmat*p;
            [xpoint, ~, ~, found] = cspice_sincpt(method, target, et,...
                targetframe, 'NONE', sc, targetframe, p_body);
            if found
                [~, lon, lat] = cspice_reclat(xpoint);
                grid_topo{i, j} = [lon*cspice_dpr, lat*cspice_dpr];
            end
        end
    end
end
end