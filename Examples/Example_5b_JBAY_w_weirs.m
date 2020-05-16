% Example_5b_JBAY_w_weirs
% Mesh the New York Jamaica bay (JBAY) region in
% high resolution with two 15-30 m wide weirs at the mouth of the estuary.
clc; clearvars

addpath(genpath('../utilities/'));
addpath(genpath('../datasets/'));
addpath(genpath('../m_map/'));

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
%% The weirs input is a structure with the following options
load weirs_struct.mat
% Weirs is an array of
%  structs each with fields:
%              X: [N—1 double] % x georgraphic coordinates of crestline
%              Y: [N—1 double] % y geographic coordinates of crestline
%          width: 5 % seperation of front and back face in meters
%        min_ele: 20 % minimum resolution along faces of weir in meters
%    crestheight: 5 % in meters
% weirs(1).X = [
%   -73.9396
%   -73.9304
%   -73.9289];
% weirs(1).Y = [
%    40.5756
%    40.5548
%    40.5507];
% weirs(1).width = 10;
% weirs(1).min_ele = 100;
% weirs(1).crestheight = 5;


% NOTE: weirs are modeled as "inner" geometry or artifical islands. A thin
% knife edge is added to both tips of the weirs to avoid bad element transistions.
% to nearby non-weir elements.
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
m = make_bc(m,'auto',gdat,'both',[],5); % auto make_bc with depth_lim = 5
m.bd = [];                              % delete land bcs and replace with weirs below
m = make_bc(m,'weirs',gdat,1,[weirs.crestheight]);  % add bcs for the weirs

%% STEP 6: interpolate bathy and plot and save the mesh
m = interp(m,gdat,'nan','fill','mindepth',1); % interpolate bathy to the
                % mesh with fill nan option to make sure corners get values
plot(m,'bd'); % visualize your boundaries
plot(m,'bmesh'); % plot triangulation and bathy
caxis([-10 0]) ;
save('JBAY_HR.mat','m')
