% Example_12_Riverineflow: Guiding for make a riverine flow input file f.20

clearvars; clc;

addpath(genpath('utilities/'));
addpath(genpath('datasets/'));
addpath(genpath('m_map/'));

%% STEP 1: set mesh extents and set parameters for mesh. 
%% The greater US East Coast and Gulf of Mexico region

bbox = [112 20; 112 24; 116 24; 116 20; 112 20]; %polygon boubox
min_el    = 1e3;  		        % minimum resolution in meters.
max_el    = 20e3; 		        % maximum resolution in meters. 
dt        = 0;                  % Automatically set timestep based on nearshore res
grade     = 0.3;               % mesh grade in decimal percent. 
R         = 3; 			        % Number of elements to resolve feature.
  
%% STEP 2: specify geographical datasets and process the geographical data
%% to be used later with other OceanMesh classes...
dem       = 'SRTM15+V2.1.nc';
coastline = 'GSHHS_f_L1';
gdat = geodata('shp',coastline,'dem',dem,'h0',min_el,'bbox',bbox);
plot(gdat,'shp'); 
plot(gdat,'dem');
            
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,'fs',R,'max_el',max_el,'dt',dt,'g',grade);
         
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','lam');
mshopts = mshopts.build; 

%% Plot and save the msh class object/write to fort.14
m = mshopts.grd; % get out the msh object
m = interp(m,gdat,'mindepth',1); % interpolate bathy to the mesh with minimum depth of 1 m
m = lim_bathy_slope(m,0.1,0);
m = make_bc(m,'auto',gdat);               % make the nodestring boundary conditions

dt = fh.dt;
CFL = CalcCFL(m,dt); 
max(CFL)
min(CFL)
m = BoundCr(m,dt);

plot(m,'type','bd');  % plot mesh on native projection with boundary conditions
plot(m,'type','b');   % plot bathy on native projection
save('Riverine.mat','m'); write(m,'Riverine');