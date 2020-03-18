% Example_5b_JBAY_w_weirs
% Mesh the New York Jamaica bay (JBAY) region in 
% high resolution with two 15-30 m wide weirs at the mouth of the estuary. 
clc; clearvars

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
load weirs_struct
gdat = geodata('shp',coastline,...
               'dem',dem,...
               'bbox',bbox,...
               'h0',min_el,...  
               'weirs',weirs);
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,...
            'fs',R,...
            'dt',dt,...
            'max_el',max_el,...
            'g',grade);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','trans');
% now build the mesh with your options and the edge function.
mshopts = mshopts.build; 

m = mshopts.grd;

%% STEP 5: Manually specify open boundaries and weir crest heights
m = makens(m,'auto',gdat,[],5);       % auto makens
m.bd = [];                            % reset land boundaries
m = makens(m,'weirs',gdat,1,[weirs.crestheight]);  % add nodestrings for the weirs

%% STEP 6: interpolate bathy and plot and save the mesh
m = interp(m,gdat,'nan','fill','mindepth',1); % interpolate bathy to the 
                % mesh with fill nan option to make sure corners get values
plot(m,'bd'); % visualize your boundaries 
plot(m,'bmesh'); % plot triangulation and bathy
caxis([-10 0]) ; 
save('JBAY_HR.mat','m')
