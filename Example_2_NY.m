% Example_2_NY: Mesh the New York region in high resolution
addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))
%% STEP 1: Set mesh extents and set parameters for mesh.
bbox = [-74.5 -73.8 	% lon_min lon_max
        40.5 40.9]; 	% lat_min lat_max
min_el    = 30;  	% Minimum resolution in meters.
max_el    = 1e3; 	% Maximum resolution in meters. 
max_el_ns = 240;        % Maximum resolution nearshore in meters.
dt        = 2;          % Encourage mesh to be stable at a 2 s timestep
grade     = 0.20; 	% Mesh grade in decimal percent.
R         = 3;    	% Number of elements to resolve feature width.
%% STEP 2: Specify geographical datasets and process the geographical data 
%% to be used later with other OceanMesh classes.
coastline = 'PostSandyNCEI';
dem       = 'PostSandyNCEI.nc';
gdat = geodata('shp',coastline,...
               'dem',dem,...
               'bbox',bbox,...
               'h0',min_el); 
%% STEP 3: Create an edge function class.
fh = edgefx('geodata',gdat,...
            'fs',R,'max_el_ns',max_el_ns,...
            'max_el',max_el,'dt',dt,'g',grade);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',{fh},'bou',{gdat},'nscreen',1,'plot_on',1,'itmax',150);
mshopts = mshopts.build; 
%% STEP 5: Plot it and write a triangulation fort.14 compliant file to disk.
plot(mshopts.grd,'tri');
write(mshopts.grd,'NY_HR');
