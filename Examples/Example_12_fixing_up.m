% Example_12_remeshing
% 
% This example highlights a workflow on how to re-mesh and insert patches
% using user defined mesh sizing functions from OceanMesh2D preserving
% subdomain boundaries exactly.
%
% By Keith Roberts, 2020, USP

clearvars; clc;

addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../m_map/'))

%% STEP 1: Here we mimic the steps that occurred in Example_1_NZ.m
bbox = [166 176;		% lon_min lon_max
    -48 -40]; 		% lat_min lat_max
min_el    = 1e3;  		% minimum resolution in meters.
max_el    = 100e3; 		% maximum resolution in meters.
max_el_ns = 5e3;        % maximum resolution nearshore in meters.
grade     = 0.35; 		% mesh grade in decimal percent.
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
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'nscreen',5,'proj','trans');
mshopts = mshopts.build;

%% STEP 5: Plot it and write a triangulation fort.14 compliant file to disk.
% Get out the msh class and put on nodestrings
m = mshopts.grd;

%% Extract the region which you want to remesh
% For the purpose of this example, we have drawn a random polygon on top
% of the mesh for which we would like to remesh.
hole = [  174.1578  -45.9624
    172.5435  -46.8280
    171.7715  -46.1261
    172.5903  -45.0266
    174.4151  -44.5821];
% This extracts the parent mesh but with the polygon removed.
subdomain = extract_subdomain(m,hole,1);
% Visualize the subdomain.
plot(subdomain);
%% Get the boundary of this hole in the mesh
% Follow the instructions in the title of the plot and click on one of the 
% points on the hole's boundary and type the number into the screen.
poly = get_poly(subdomain);
%% Re-meshing 
% Using the polygon we just extracted (poly) and the original mesh size 
% function (fh), remesh this hole with Delaunay refinement mesh2d. 
% Recall poly is a cell-array so pass it the hole you want to mesh!

subdomain_new = mesh2dgen(poly{1},fh); 

%% Merge the patch back into the parent mesh.
% Since the boundaries match identically, this is a trivial merge
% and we use the `match` option.
m_new = plus(subdomain_new, m, 'match');

% Visualize the final result (the red indicates the area that we remeshed).
plot(m_new); hold on; plot(poly{1}(1:end-1,1),poly{1}(1:end-1,2),'r-','linewi',3);


