% Example_10_Multiscale_Smoother:
% An idealized test for multiscale nesting using boxes with a large min_el
% ratio
clearvars; clc; close all

PREFIX = '10_Multiscale_Smoother';

%% STEP 1: set mesh extents and set parameters for mesh.
bbox = [0 1; 0 1];
boubox = [0 0; 1 0; 1 1; 0 1; 0 0; NaN NaN];
min_el = 1e3;

bbox2 = [-1 2; -1 2];
boubox2 = [-1 -1; 2 -1; 2 2; -1 2; -1 -1; NaN NaN];
min_el2 = min_el*10;

%% STEP 2: specify geographical datasets and process the geographical data
% to be used later with other OceanMesh classes...
gdat1 = geodata('pslg',boubox,'bbox',bbox,'h0',min_el);
gdat2 = geodata('pslg',boubox2,'bbox',bbox2,'h0',min_el2);

%% STEP 3: create an edge function class
fh1 = edgefx('geodata',gdat1,'g',0.2);
fh2 = edgefx('geodata',gdat2,'g',0.2);

%% STEP 4: Pass your edgefx class object along with some meshing options
% and build the mesh...
mshopts = meshgen(...
    'ef',{fh2, fh1},...
    'bou',{gdat2,gdat1},...
    'plot_on',1,...
    'qual_tol',0.0025,...
    'cleanup',0);
mshopts = mshopts.build;

m = mshopts.grd;

% Plotting mesh without cleaning
plot(m)
title('Mesh without cleaning')

% clean
m = m.clean;

% Plotting mesh with cleaning
plot(m)
title('Mesh after cleaning')

%% Export plots
figs = get(0,'children');
for f = 1:numel(figs)
    fname = sprintf('%s_plot%i',PREFIX,figs(f).Number);
    print(figs(f).Number,fname,'-dpng','-r200');
end

%% Save mesh files
% Save as a msh class
save(sprintf('%s_msh.mat',PREFIX),'m');
% Write an ADCIRC fort.14 compliant file to disk.
write(m,sprintf('%s_mesh',PREFIX))
