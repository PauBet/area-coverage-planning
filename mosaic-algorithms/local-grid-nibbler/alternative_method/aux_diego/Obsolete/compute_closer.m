% Compute area closer corner to the zero target
function closer_corner = compute_closer(areapoints,zerotarget)

    % last point in areapoints is the same as the first (for a better plotting)
    diff_vecs = areapoints(1:end-1,:)-zerotarget.';
    diff = vecnorm(diff_vecs.');
    [~, closer_index] = min(diff);
    closer_corner = areapoints(closer_index,:);

end