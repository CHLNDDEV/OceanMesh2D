% Example_6c_GBAY_w_floodplain2
% Continue on from Example_6_GBAY.m by building on
% a floodplain onto the mesh using Mesh2D front propagation technique
% NOTE: need to download mesh2d separately
% (https://github.com/dengwirda/mesh2d)
%%
clearvars; close all; clc;

addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))

%% STEP 1: set mesh extents and set parameters for mesh.
min_el    = 60;  		    % Minimum mesh resolution in meters.
max_el    = [1e3 0 -inf     % Globally, maximum mesh resolution in meters.
             250 inf 0];    % Overland, maximum mesh resolution in meters.
grade     = [0.25 0 -inf    % Use a spatially variable gradation rate overland.
             0.05 inf 0] ;
angleOfReslope = 60 ;       % Control width of channel by changing angle of reslope.
ch = 0.1 ;                  % Scale resolution propotional to depth nearby thalweg.
fs = 3 ;                    % Place 3 vertices per width of shoreline feature.

%% STEP 2: specify geographical datasets and process the geographical data
%% to be used later with other OceanMesh classes...
coastline = 'us_medium_shoreline_polygon';
demfile   = 'galveston_13_mhw_2007.nc';

gdatuw = geodata('shp',coastline,...
                 'dem',demfile,...
                 'h0',min_el);

load ECGC_Thalwegs.mat % Load the Channel thalweg data

%% STEP 3: create an edge function class
fh = edgefx('geodata',gdatuw,...
            'fs',fs,...
            'ch',ch,...
            'AngOfRe',angleOfReslope,...% control the width
            'Channels',pts2,...
            'g',grade,...
            'max_el',max_el);  
        
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',fh,'bou',gdatuw,'plot_on',1,'proj','lambert');
% now build the mesh with your options and the edge function.
mshopts = mshopts.build;
muw = mshopts.grd ;

%% STEP 5: Specify the outer PSLG to mesh up to and stitch to the ocean mesh polygon
gdatfp.outer = gdatuw.boubox(1:end-2,:);
% need to make sure that we start and end near shoreline-bbox edge
gdatfp.outer(1,2) = 28.979; gdatfp.outer(end,2) = 29.579; 
% the ocean side intersects with bbox, not ideal.. but we need to expand
% the floodplain out sligthly.. 
gdatfp.outer(1:2,1) = gdatfp.outer(1:2,1) - max_el(2,1)/111e3;
[la,lo] = interpm(gdatfp.outer(:,2),gdatfp.outer(:,1),max_el(2,1)/111e3);
gdatfp.outer = [lo, la];
polygon = Stitch_Shoreline_to_Upper_Contour(muw,gdatfp);

%% STEP 6: Make the floodplain mesh with Mesh2D
mfp = mesh2dgen(polygon,fh);

%% STEP 7: Interpolate bathy and merge
% limits ocean side to at least 1 m depth
muw = interp(muw,demfile,'mindepth',1);
% uses wider cell-average on floodplain, and ensures at least 1 ft above ground 
mfp = interp(mfp,demfile,'N',3,'maxdepth',-0.3,'nan','fill');
[~,mfp] = setProj(mfp,1,'lambert',1);
m = plus(muw,mfp,'match');

%% STEP 8: Plot bathy and mesh to check
plot(m,'bmesh')
demcmap([-20 10]);