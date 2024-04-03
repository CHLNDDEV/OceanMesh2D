% Example_8_AK: Mesh the western Alaska and Bering Shelf region 
% 
% Note about the example:
% The domain crosses the -180/180 meridian but does not wrap
% completely around so it requires the user to specify the bbox in 0-360
% deg format to avoid the confusion. The geodata routine handles this 
% internally even though the shapefile (and/or the dem) is in 
% -180 to 180 deg format. 

clearvars; clc; close all

PREFIX = '8_AK';

%% STEP 1: set mesh extents and set parameters for mesh.
% This is a polygonal shape for west Alaska/Bering Shelf 
% in 0-360 deg format
load('ak_outerpoly.mat');

min_el = 5e3;       % minimum resolution in meters.
max_el = 50e3;      % maximum resolution in meters. 
grade = 0.25;       % mesh grade in decimal percent.

%% STEP 2: specify geographical datasets and process the geographical data 
% to be used later with other OceanMesh classes...
dem = 'SRTM15+.nc';
coastline = 'GSHHS_f_L1';

gdat = geodata(...
    'shp',coastline,...
    'bbox',bbox,...
    'h0',min_el);

% plotting the gdat to show it crosses the 180/-180 and 
% is in 0-360 format meridian. 
plot(gdat,[],'lam')

%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,'dis',grade,'max_el',max_el,'g',grade);
        
%% STEP 4: Pass your edgefx class object along with some meshing options 
% and build the mesh...
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'nscreen',5,'proj','lam');
mshopts = mshopts.build; 

%% STEP 5: Extract the msh class, make the boundary conditions and plot
m = mshopts.grd; % Get out the msh class
m = make_bc(m,'auto',gdat,'distance'); % make the boundary conditions
plot(m,'type','bd'); % plot the boundary conditions on the mesh

m = interp(m,dem,'mindepth',5);

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
