% Example_9_edgefx_filter.m
% 
%   This script serves as an example of how the edgefx_filter function can
% ensure that the smallest possible target edgelength function in a series
% of edgefunctions is carried through to the innermost edgefunction before
% being used to generate a mesh in meshgen. The reason that this function
% can sometimes be necessary is because in different locations, different
% controlling edgefunction parameters can result in larger target
% edgelengths being put "above" or "inside" the edgefunction that meshgen
% utilizes. A simple example, and one used here, is if the channels
% parameter is utilized in an outer edgefunction but not an inner
% edgefunction. As can be seen by the two meshes generated using this
% script, the edgefx_filter fixes this problem and ensures the channels are
% still incorporated in the final mesh.
%
% Coleman Blakely
% September 4, 2020
%% Setup
clear;clc;hold off;close all;
%% Add libraries
addpath(genpath('C:\Users\Diane\Documents\GitHub\OceanMesh2D\'))
%% Define data sources
% Shoreline
shoreline = 'GSHHS_f_L1';
% Bathy data
dem = 'SRTM15+V2.nc';
%% Define meshing domain and parameters
bbox{1} = [-83, -82; 27 28.5];        % outermost bbox
min_el = 50;                       % minimum edgelength
max_el = 1e3;                       % maximum edgelength
grade = 0.35;                       % maximum change in size of edgelength
R = 3;                              % elements required to resolve a feature
load ECGC_Thalwegs.mat
%% Build geodata and edgelength function for outermost domain
gdat{1} = geodata('shp',shoreline,'dem',dem,'bbox',bbox{1},'h0',min_el);
fh{1} = edgefx('geodata',gdat{1},'fs',R,'max_el',max_el,'g',grade,'Ch',1,'Channels',pts2);
%% Build inner domain
% To show how the edgefunction filter works, build another edgefunction
% within the outer domain with a larger minimum edgelength
% min_el = 150;
bbox{2} = [-82.8, -82.4; 27.25 28.25];
gdat{2} = geodata('shp',shoreline,'dem',dem,'bbox',bbox{2},'h0',min_el);
fh{2} = edgefx('geodata',gdat{2},'fs',R,'max_el',max_el,'g',grade);
% plotedgefx(fh,1);
% bound = bbox2bound(bbox{2});
% plot(bound(:,1),bound(:,2),'r-')
%% Generate the mesh using the unfiltered edgefunctions
m_unfiltered = meshgen('ef',fh,'bou',gdat,'plot_on',0);
m_unfiltered = m_unfiltered.build;
m_unfiltered = m_unfiltered.grd;
%% Use the edgefunction filter
fh = edgefx_filter(fh,gdat);
% plotedgefx(fh,2)
%% Generate the mesh using the filtered edgefunction
m_filtered = meshgen('ef',fh,'bou',gdat,'plot_on',0);
m_filtered = m_filtered.build;
m_filtered = m_filtered.grd;
%% Plot the two meshes
plot(m_unfiltered,'tri',0)
title('Unfiltered Edgefunction')
plot(m_filtered,'tri',0)
title('Filtered Edgefunction')