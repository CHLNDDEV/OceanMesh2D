% convert Admesh (ADCIRC) *.14 in Cartesin coordinates to SMS *.2dm using Matlab on Windows
% Author: Jun Sasaki  Coded on September 20, 2021  Updated on September 20, 2021

clearvars; clc;

% Set the following paths
if ispc % Code to run on Windows
    top_dir = 'D:/Github/'
elseif isunix % Code to run on Linux
    top_dir = '~/Github/'
end

% Specify input Admesh 14 file
base_dir = [top_dir, '/OceanMesh2D/fvcom/']  % Specify the directory containing admesh_file
admesh_file = 'dem_00_01_change' % without extension

addpath([top_dir, '/fvcom-toolbox/utilities/']);
addpath([top_dir, '/fvcom-toolbox/fvcom_prepro/']);

Mobj=read_admesh_mesh('msh',fullfile(base_dir,[admesh_file,'.14']),'coordinate','cartesian');
Mobj

write_SMS_2dm(fullfile(base_dir,[admesh_file,'.2dm']),Mobj.tri,Mobj.x,Mobj.y,Mobj.h,[])
