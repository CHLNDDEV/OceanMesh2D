%Example_2_NY: Mesh the New York Jamaica bay (JBAY region in high resolution.
addpath(genpath('utilities/'));
addpath(genpath('datasets/'));
addpath(genpath('m_map/')); 
%% STEP 1: set mesh extents and set parameters for mesh.
bbox      = [-73.97 -73.75 	% lon_min lon_max
              40.5 40.68]; 	% lat_min lat_max
min_el    = 15.0;  		    % Minimum resolution in meters.
max_el    = 1e3; 
dt        = 2;              % Ensure mesh is stable at a 2 s timestep
grade     = 0.15;           % Mesh grade in decimal percent.
R         = 3;              % Number of elements to resolve feature width.
%% STEP 2: specify geographical datasets and process the geographical data 
%% to be used later with other OceanMesh classes...
coastline = 'PostSandyNCEI'; 
dem       = 'PostSandyNCEI.nc'; 
gdat = geodata('shp',coastline,...
               'dem',dem,...
               'bbox',bbox,...
               'h0',min_el);
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,...
            'fs',3,...
            'dt',dt,...
            'max_el',max_el,...
            'g',grade);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',{fh},'bou',{gdat},'plot_on',1);
% now build the mesh with your options and the edge function.
mshopts = mshopts.build; 
%% STEP 5: Plot it 
plot(mshopts.grd,'tri');
