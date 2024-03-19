function loadKernels(metakr)

% Clean pool of kernels
endSPICE

% Get homespice
HOMESPICE = getHomeSpice;

% Load kernels
for i=1:length(metakr)
    if mod(i,2)==1
        url = metakr{i};
        continue;
    end
    kfile = fullfile(HOMESPICE, 'kernels', metakr{i});
    if exist(kfile,'file')~=2
        websave(kfile,url);
    end
    cspice_furnsh(kfile);
end

% Total loaded kernels
fprintf("Kernel pool: %d\n", cspice_ktotal('all'));

end