% Example_9_TBay.m (Tampa Bay)
% 
% This script serves as an example of how the enforce_min_ef function can
% ensure that the smallest possible target edgelength function in a series
% of edgefunctions is carried through to the innermost edgefunction before
% being used to generate a mesh in meshgen. The reason that this function
% can sometimes be necessary is because in different locations, different
% controlling edgefunction parameters can result in larger target
% edgelengths being put "above" or "inside" the edgefunction that meshgen
% utilizes. A simple example, and one used here, is if the channels
% parameter is utilized in an outer edgefunction but not an inner
% edgefunction. As can be seen by the two meshes generated using this
% script, the enforceMin option fixes this problem and ensures the channels
% are still incorporated in the final mesh.
%
% Limitation: The min_el (h0) in the inner cannot be larger than in the outer
%
% Coleman Blakely
% September 4, 2020
% Modifications by William pringle Sep 19, 2020

%% Setup
clearvars; clc; close all;

%% Define data sources
% Shoreline
shoreline = 'GSHHS_f_L1';
% Bathy data
dem = 'SRTM15+V2.nc';
%% Define meshing domain and parameters
bbox{1} = [-83, -82; 27 28.5];      % outermost bbox
min_el = 100;                       % minimum edgelength
max_el = 1e3;                       % maximum edgelength
grade = 0.35;                       % maximum change in size of edgelength
R = 3;                              % elements required to resolve a feature
load ECGC_Thalwegs.mat

%% Build geodata and edgelength function for outermost domain
gdat{1} = geodata('shp',shoreline,'dem',dem,'bbox',bbox{1},'h0',min_el);
fh{1}   = edgefx('geodata',gdat{1},'fs',R,'max_el',max_el,...
                 'g',grade,'Ch',1,'Channels',pts2);

%% Build inner domain
% To show how the enforce_min_ef function works, build another edgefunction
% within the outer domain that does not use the channel function
bbox{2} = [-82.8, -82.4; 27.25 28.25];
gdat{2} = geodata('shp',shoreline,'dem',dem,'bbox',bbox{2},'h0',min_el);
fh{2}   = edgefx('geodata',gdat{2},'fs',R,'max_el',max_el,'g',grade);

%% Generate the mesh without enforcing the min of all edgefunctions
m_nomin = meshgen('ef',fh,'bou',gdat,'enforceMin',0,'plot_on',0);
m_nomin = m_nomin.build;
m_nomin = m_nomin.grd;

%% Generate the mesh with enforcing the min of all edgefunctions 
%-> by default enforceMin = 1 so do not need to set usually
m_min = meshgen('ef',fh,'bou',gdat,'enforceMin',1,'plot_on',0);
m_min = m_min.build;
m_min = m_min.grd;

%% Plot the two meshes
plot(m_nomin,'tri')
title('Without Enforcing Minimum of all Edgefunctions')
plot(m_min,'tri')
title('Enforcing Minimum of all Edgefunctions')