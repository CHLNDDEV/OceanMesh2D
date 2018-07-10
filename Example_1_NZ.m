% Example_1_NZ: Mesh the South Island of New Zealand
addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))
%% STEP 1: set mesh extents and set parameters for mesh.
bbox = [166 176 		% lon_min lon_max
        -48 -40]; 		% lat_min lat_max
min_el    = 250;  		% minimum resolution in meters.
max_el    = 20e3; 		% maximum resolution in meters. 
max_el_ns = 1e3;        % maximum resolution nearshore in meters.
grade     = 0.20; 		% mesh grade in decimal percent.
R         = 3;    		% number of elements to resolve feature width.
%% STEP 2: specify geographical datasets and process the geographical data 
%% to be used later with other OceanMesh classes...
coastline = 'GSHHS_f_L1';
gdat = geodata('shp',coastline,'bbox',bbox,'h0',min_el);
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,...
            'fs',R,'max_el_ns',max_el_ns,...
            'max_el',max_el,'g',grade);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopt = meshgen('ef',{fh},'bou',{gdat},'nscreen',5,'plot_on',1,'itmax',75);
mshopt = mshopt.build; 
%% STEP 5: Plot it and write a triangulation fort.14 compliant file to disk.
plot(mshopt.grd,'tri');
write(mshopt.grd,'South_Island_NZ');
