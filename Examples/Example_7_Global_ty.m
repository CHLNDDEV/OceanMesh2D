% Example_7_Global: Make a global mesh
clearvars; clc;

addpath('..')
addpath(genpath('../utilities/'));
addpath(genpath('../datasets/'));
addpath(genpath('../m_map/'));

%% STEP 1: set mesh extents and set parameters for mesh. 
%% The greater US East Coast and Gulf of Mexico region
%min_el    = 4e3;  	             % minimum resolution in meters.
min_el    = 2.25e3;  	             % minimum resolution in meters.
bbox      = [-180 180; -88.5 90]; % lon min lon max; lat min lat max
%max_el    = 20e3; 		         % maximum resolution in meters. 
max_el    = 25e3; 		         % maximum resolution in meters. 
wl        = 30;                  % 30 elements resolve M2 wavelength.
dt        = 0;                   % Only reduces res away from coast
%grade     = 0.25;                % mesh grade in decimal percent. 
grade     = 0.35;                % mesh grade in decimal percent. 
%R         = 3; 			     % Number of elements to resolve feature.
%slp       = 10;                  % slope of 10
slp       = 10;                  % slope of 10
fl        = -20;   % added 
dis       = grade; % added

%outname = 'Global_4km_20km';
outname = 'Global_2.2km_25km';

%load('global_gdat_fh')

% %% STEP 2: specify geographical datasets and process the geographical data
% %% to be used later with other OceanMesh classes...
dem       = 'SRTM15+.nc';
coastline = {'GSHHS_f_L1','GSHHS_f_L6'};
gdat1 = geodata('shp',coastline,'dem',dem,...
                'bbox',bbox,'h0',min_el);
%             
% %% STEP 3: create an edge function class
% %fh1 = edgefx('geodata',gdat1,...
% %             'fs',R,'wl',wl,'max_el',max_el,...
% %             'slp',slp,'dt',dt,'g',grade); 
fh1 = edgefx('geodata',gdat1,...
             'dis',dis,'wl',wl,'max_el',max_el,'slp',slp,...
             'g',grade); 
% % 
% save('global_gdat_fh','-v7.3','gdat1','fh1');
               
%% STEP 4: Pass your edgefx class object along with some meshing options 
%% and build the mesh...
mshopts = meshgen('ef',fh1,'bou',gdat1,'nscreen',10,'plot_on',1,...
                  'proj','stereo','dj_cutoff',3.5e5);
mshopts = mshopts.build; 

%% STEP 5: Match points and edges across boundary
m = mshopts.grd; % get out the msh object 
% Plotting the triangulation on Robinson projection
plot(m,'type','tri','proj','Robinson');
print('Robinson_mesh_tri','-dpng','-r300')

save([outname '.mat'],'m'); 

%write(m,outname);
