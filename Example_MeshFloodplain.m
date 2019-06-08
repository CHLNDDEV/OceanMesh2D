%%
% Example_MeshFloodplain: Continue on from Example_6_GBAY.m by building on
% a floodplain onto the mesh.
%%
clearvars; close all; clc;

addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))
%% STEP 1: set mesh extents and set parameters for mesh.
bbox = [-95.40 -94.4;
         29.14  30.09];
min_el    = 60;  		    % Minimum mesh resolution in meters.
max_el    = [1e3 0 -inf     % Globally, maximum mesh resolution in meters.
             1e3 inf 0];    % Overland, maximum mesh resolution in meters.
grade     = [0.25 0 -inf    % Use a spatially variable gradation rate overland.
             0.05 inf 0] ;
angleOfReslope = 60 ;       % Control width of channel by changing angle of reslope.
ch = 0.1 ;                  % Scale resolution propotional to depth nearby thalweg.
fs = 3 ;                    % Place 3 vertices per width of shoreline feature. 
%% STEP 2: specify geographical datasets and process the geographical data
%% to be used later with other OceanMesh classes...
coastline = 'us_medium_shoreline_poly';
demfile   = 'galveston_13_mhw_2007.nc';

gdat = geodata('shp',coastline,...
    'dem',demfile,...
    'bbox',bbox,...
    'h0',min_el);

load ECGC_Thalwegs.mat % Load the Channel thalweg data
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,...
    'fs',fs,...
    'ch',ch,...
    'AngOfRe',angleOfReslope,...% control the width
    'Channels',pts2,...
    'g',grade,...
    'max_el',max_el);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','lambert');
% now build the mesh with your options and the edge function.
mshopts = mshopts.build;
%% STEP 5: Get fixed constraints and update gdat with overland meshing domain.
muw = mshopts.grd ;
muw = makens(muw,'auto',gdat) ; % apply so that extractFixedConstraints only grabs the shoreline constraints.

[pfix,egfix] = extractFixedConstraints(muw) ;

% 10-m contour extracted from the Coastal Relief Model.
coastline = 'us_coastalreliefmodel_10mLMSL';
demfile   = 'galveston_13_mhw_2007.nc';

gdat = geodata('shp',coastline,...
    'dem',demfile,...
    'bbox',bbox,...
    'h0',min_el);

gdat.inpoly_flip =  mod(1,gdat.inpoly_flip) ; % if the meshing domain is inverted, you can always flip it .

% Here we pass our constraints to the mesh generator (pfix and egfix). 
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','lambert','pfix',pfix,'egfix',egfix);

% now build the mesh with your options and the edge function.
mshopts = mshopts.build;

m = mshopts.grd ;


