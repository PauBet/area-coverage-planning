clc
clear 
close all

tic

% SPICE KERNELS LOAD 

for a = 1:1
    addpath(genpath('E:\RESSlib'));
    % A LIBRARY CALLED CLIPPER2 HAS BEEN USED IN ORDER TO MAKE FASTEST GEOMETRY
    % COMPUTATIONS SUCH AS INTERSECTIONS OF POLYGONS.
    addpath(genpath('E:\MATLABR2022\toolbox\clipper2'));   
    METAKR={ 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls', ...
            'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/gll_951120_021126_raj2007.bsp', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/s030916a.bsp', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/s000131a.bsp', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/s980326a.bsp', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/spk/s970311a.bsp', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/pck/pk96030a.tpc', ...
            'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/kernels/sclk/mk00062a.tsc', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/misc/fk/gll_v0.tf', ...
            'https://naif.jpl.nasa.gov/pub/naif/GLL/misc/ibelgacem/gll_plt_inst_ikfk_v01.tf', ...
            'pck00010_msgr_v23_europa2020.tpc',...
            'cke15ahb_plt.bc',...
            'cke15ahb_plt.bc',... % E15
            'cke17f_plt.bc',...
            'cke17bfb_plt.bc',... % E17
            'cke17afc_plt.bc',...
            'cke19ahc_plt.bc',... % E19
            'cke19f_plt.bc',...
            'cke19cdc_plt.bc',...
            'cke19bdh_plt.bc',...
            'cke14f_plt.bc',...  % E14
            'cke14aje_plt.bc',...
            'cke12alc_plt.bc',... % E12
            'cke12f_plt.bc',...
            'ckg07pbq_plt.bc',... % G7
            'ckg07rtq_plt.bc',...
            'ckg07b_plt.bc',...
            'ckg07a1g_plt.bc',...
            'ckg07a2h_plt.bc',...
            'cke06ape_plt.bc',... % E6
            'cke06rtq_plt.bc',...
            'cke06pbq_plt.bc',...
            'cke06b_plt.bc',...
            'cke06cge_plt.bc',...
            'cke04pbq_plt.bc',... % E4
            'cke04rtq_plt.bc',...
            'cke04b_plt.bc',...
            'cke04ckj_plt.bc',...
            'cke04bii_plt.bc',...
            'cki25bde_plt.bc',... %I25
            'cki25adf_plt.bc',...
            'cki25f_plt.bc',...
            'cke26f_plt.bc',... %E26
            'cke26afc_plt.bc',... 
            'ckee11v3_plt.bc',... % E11
            'cke11aig_plt.bc',...
            'cke11pbq_plt.bc',...
            'cke11b_plt.bc',...
            'cke11rtq_plt.bc',...
            'ckc03pbq_plt.bc',... %C3
            'ckc03rtq_plt.bc',... 
            'ckc03clb_plt.bc',...
            'ckc03ah1_plt.bc',...
            'ckc03ak2_plt.bc',...
            'ckc03bhb_plt.bc',...
            'ckc03b_plt.bc'};
            %'galssi_eur_usgs2020.bc'};       
    initSPICEv(fullK(METAKR));
end

%% DATA AND PARAMETERS

% INSTANTS TO TEST
% E6:  1997-02-20T16:08:42.514 --> ok
% E11: 1997-11-06T19:24:04.407 --> ok
% E12: 1997-12-16T11:11:36.136 --> ok
% E14: 1998-03-29T13:17:12.815 --> ok
% E15: 1998-05-31T20:10:11.556 --> ok
%      1998-05-31T21:07:15.955 --> ok
%      1998-05-31T22:15:34.218 --> ok
% E17: 1998-09-26T02:50:38.475 --> ok
%      1998-09-26T03:41:55.137 --> ok
%      1998-09-26T04:56:27.138 --> ok
% E19: 1999-02-01T01:17:16.465 --> ok

et0 =  cspice_str2et('1998-05-31T20:10:11.556'); %Starting time for the mosaic
target_body = 'EUROPA'; %Target name
abcCorr = 'NONE'; %Aberration correction
sc = '-77'; %Observer code
inst = 'GLL_SSI'; %Instrument code
target_fixed = append('IAU_',target_body); %Target reference frame
delta_T = 2; %Time step between obersvations
Body_rad = cspice_bodvrd(target_body,'RADII',3); % BODY RADIUS

%% PARAMETERS USED FOR THE cspice_limbpt FUNCTION
ref_Vec = [ 0.0, 0.0, 1.0 ].';
NMETH  = 2;
MAXN   = 100;
schstp = 1.0d-4;
soltol = 1.0d-7;
ncuts  = 150;
delrol = cspice_twopi() / ncuts;

%% PLOT THE LIMB
[~, limb, ~, ~] = cspice_limbpt('TANGENT/ELLIPSOID',target_body,et0,target_fixed,abcCorr,'CENTER',sc,ref_Vec,delrol,ncuts,schstp,soltol,ncuts);
[~, limb_lon, limb_lat] = cspice_reclat(limb);
limb_lon = limb_lon*(180/pi);
limb_lon(limb_lon<0) = limb_lon(limb_lon<0)+360;
limb_lat = limb_lat*(180/pi);
plot([limb_lon, limb_lon(1)],[limb_lat, limb_lat(1)])

%% FIGURE PREPROCESSING (TO CHOOSE THE ROI) %%%%%%%%%%%
hold on
xlim([0 360]);
ylim([-90 90]);
xlabel('longitude [°]','FontSize', 12)
ylabel('latitude [°]','FontSize', 12)

%% GENERATE FOOTPRINT FUNCTION AND PLOT FIRST FOOTPRINT

% Zero Target
steps_zero = 40; % Number of points for the footprint function

%FOOTPRINT FUNCTION
footprint_func = @(x,y,z) footprint_coverage_final(x, sc, inst, y, target_body, target_fixed,z); 

%ZERO CONDITIONS
[zerofootprint,~, ] = footprint_func(et0,0,steps_zero);

toc

zerotarget = [mean(zerofootprint(1,:));
              mean(zerofootprint(2,:))];
plot(zerotarget(1),zerotarget(2),'.', 'MarkerSize', 9, 'Color', 'red')
hold on

%% INTRODUCING ROI MANUALLY
% add area of interest 
n = 10; % number of points to define the area
[areapoints_real] = add_area_manual(n);

tic

tobs = 2;

[real_targets] = neighbour_placement(et0, tobs, inst, sc, ...
    target_body, areapoints_real);

toc

% Plot results 

plot_area_and_footprints(et0,real_targets,areapoints_real,inst,target_fixed,footprint_func)


endSPICE;
