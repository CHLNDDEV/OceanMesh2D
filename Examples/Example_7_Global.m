% Example_7_Global: Make a global mesh 4km_20km

clearvars; clc; close all

PREFIX = '7_Global_4km_20km';

%% STEP 1: set mesh extents and set parameters for mesh.
% The greater US East Coast and Gulf of Mexico region
bbox = [
    -180 180        % lon min lon max
    -89  90         % lat min lat max
    ];
min_el = 4e3;       % minimum resolution in meters.
max_el = 20e3;      % maximum resolution in meters.
wl = 30;            % 30 elements resolve M2 wavelength.
dt = 0;             % Only reduces res away from coast
grade = 0.25;       % mesh grade in decimal percent.
R = 3;              % Number of elements to resolve feature.
slp  = 10;          % slope of 10

%% STEP 2: specify geographical datasets and process the geographical data
% to be used later with other OceanMesh classes...
dem       = 'SRTM15+.nc';
coastline = {'GSHHS_f_L1', 'GSHHS_f_L6'};

gdat1 = geodata('shp',coastline,'dem',dem,...
    'bbox',bbox,'h0',min_el);

%% STEP 3: create an edge function class
fh1 = edgefx('geodata',gdat1,...
    'fs',R,'wl',wl,'max_el',max_el,...
    'slp',slp,'dt',dt,'g',grade);

%% STEP 4: Pass your edgefx class object along with some meshing options
% and build the mesh...
mshopts = meshgen('ef',fh1,'bou',gdat1,'nscreen',10,'plot_on',1,...
    'proj','stereo');
mshopts = mshopts.build;

%% STEP 5: Match points and edges across boundary
m = mshopts.grd; % get out the msh object

% Plotting the triangulation on Robinson projection
plot(m,'type','resolog','proj','Robinson');

%% Export plots
figs = get(0,'children');
for f = 1:numel(figs)
    fname = sprintf('%s_plot%i',PREFIX,figs(f).Number);
    print(figs(f).Number,fname,'-dpng','-r200');
end

%% Save mesh files
% Save as a msh class
save(sprintf('%s_msh.mat',PREFIX),'m');
% Write an ADCIRC fort.14 compliant file to disk.
write(m,sprintf('%s_mesh',PREFIX))
