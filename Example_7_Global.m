% Example_7_Global: Make a global mesh
clearvars; clc;
addpath(genpath('utilities/'));
addpath(genpath('datasets/'));
addpath(genpath('m_map/'));

%% STEP 1: set mesh extents and set parameters for mesh. 
%% The greater US East Coast and Gulf of Mexico region
min_el    = 4e3;  	             % minimum resolution in meters.
bbox      = [-180 180; -87  89]; % lon min lon max; lat min lat max
max_el    = 15e3; 		         % maximum resolution in meters. 
wl        = 30;                  % 30 elements resolve M2 wavelength.
dt        = 0;                   % Only reduces res away from coast
grade     = 0.25;                % mesh grade in decimal percent. 
R         = 3; 			         % Number of elements to resolve feature.
slp       = 10;                  % slope of 10

outname = 'Global_4km';

%% STEP 2: specify geographical datasets and process the geographical data
%% to be used later with other OceanMesh classes...
dem       = 'bathy_SSG_1_120.nc'; 
coastline = {'GSHHS_f_L1','GSHHS_f_L6'};
gdat1 = geodata('shp',coastline,'dem',dem,...
                'bbox',bbox,'h0',min_el);
            
%% STEP 3: create an edge function class
fh1 = edgefx('geodata',gdat1,...
             'fs',R,'wl',wl,'max_el',max_el,...
             'slp',slp,'dt',dt,'g',grade);
                
%% STEP 4: Pass your edgefx class object along with some meshing options 
%% and build the mesh...
mshopts = meshgen('ef',fh1,'bou',gdat1,...
                  'nscreen',5,'plot_on',1,'itmax',100);
mshopts = mshopts.build; 

%% STEP 5: Match points and edges across boundary
m = mshopts.grd; % get out the msh object 
% make sure not above bbox (not at 90 deg)
m.p(m.p(:,2) > bbox(2,2),2) = bbox(2,2);

% This makes sure we match left and right side and puts on the periodic bcs
m = makens(m,'auto',gdat1);
plot(m,'bd',0) %,'Mollweide');
save([outname '.mat'],'m'); %write(m,outname);