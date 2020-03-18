% Example_2_NY: Mesh the New York region in high resolution
clearvars; clc;
addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))
%% STEP 1: Set mesh extents and set parameters for mesh.
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
               'h0',min_el); 
%% STEP 3: Create an edge function class.
fh = edgefx('geodata',gdat,...
            'fs',R,'max_el_ns',max_el_ns,...
            'max_el',max_el,'dt',dt,'g',grade);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',fh,'bou',gdat,'proj','utm','plot_on',1);
mshopts = mshopts.build; 
%% STEP 5: Plot it and write a triangulation fort.14 compliant file to disk.
% Get out the msh class and put on bathy and nodestrings
m = mshopts.grd;
m = interp(m,gdat,'mindepth',1); % interpolate bathy to the mesh with minimum depth of 1 m
m = makens(m,'auto',gdat,[],5);  % make the nodestring boundary conditions 
                                 % with min depth of 5 m on open boundary
plot(m,'bd'); plot(m,'blog');    % plot triangulation, and bathy on log scale
write(m,'NY_HR');                % write to ADCIRC compliant ascii file
