function getSPKcov(object)
% This function prints the coverage windows of an object within the SPK
% loaded kernels

% Check input
if ischar(object)
    objectID = cspice_bodn2c(object);
    objectname = upper(object);
else
    objectID = object;
    objectname = cspice_bodc2n(object);
end

% Get number of SPK kernels 
count = cspice_ktotal('spk');

%
fprintf('\n%s SPK coverage:\n', objectname)

% Retrieve loaded kernel files
for i=1:count
    % Retrieve kernel path
    kpath = cspice_kdata(i, 'spk');

    % Split string to get the filename
    kpathsplit = strsplit(kpath, '/');
    kfile = kpathsplit{end};
    
    % Calculate and print the spk coverage windows
    fprintf('\nKernel %s: \n', kfile);
    cov = cspice_spkcov(kpath, objectID, 100);
    for j=1:length(cov)/2
        fprintf('From: %s\t', cspice_et2utc(cov(2*j-1), 'C', 0));
        fprintf('To: %s\n', cspice_et2utc(cov(2*j), 'C', 0));
    end
end

end