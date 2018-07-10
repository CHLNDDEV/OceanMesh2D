% Example_4_PRVI: Mesh the west North Atlantic Ocean with high resolution
% around Puerto Rico and US Virgin Islands.
clc; clearvars

addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))

%% The constant parameters for all domains
wl        = 30;         % elements to resolve M2 wavelength
dt        = 0;          % use automatic timestep 
grade     = 0.25; 		% mesh grade in decimal percent.
R         = 5;    		% number of elements to resolve feature width.
slp       = 15;         % 2*pi/number of elements to resolve slope

%%  For relatively coarse resolution west North Atlantic Ocean
bbox = [-100 -53   	    % lon_min lon_max
         5    52.5];    % lat_min lat_max
min_el = 1000;  		% minimum resolution in meters.
max_el = 10e3;       	% maximum resolution in meters. 

coastline = 'GSHHS_f_L1';
dem       = 'topo15_compressed.nc';
gdat{1} = geodata('shp',coastline,...
                  'dem',dem,...
                  'bbox',bbox,...
                  'h0',min_el);
           
fh{1} = edgefx('geodata',gdat{1},...
               'fs',R,...
               'wl',wl,...
               'slp',slp,...
               'max_el',max_el,...
               'dt',dt,...
               'g',grade);

%% For High Resolution around Puerto Rico and US Virgin Islands
coastline = ["pr_1s_0m_contour","usvi_0m_contour","sj_0contour_closed"];
    
dems      = ["pr_1s.nc","usvi_1_mhw_2014.nc", "san_juan_19_prvd02_2015.nc"];
    
for ii = 1:length(dems)
    % use same parameters as coarse mesh, just change min_el
    if ii == length(dems)
        min_el    = 10;  % minimum resolution in meters.
    else
        min_el    = 30;  % minimum resolution in meters.
    end
    % bbox is taken automatically from the DEM
    gdat{ii+1} = geodata('shp',coastline{ii},...
                         'dem',dems{ii},...
                         'h0',min_el);
    fh{ii+1} = edgefx('geodata',gdat{ii+1},...
                    'fs',R,...
                    'wl',wl,...
                    'slp',slp,...
                    'max_el',max_el,...
                    'dt',dt,...
                    'g',grade);
end

%% Pass your edgefx class objects along with some meshing options 
%% and build the mesh... 
% (note that the nested edgefxs will be smoothed together with this call)
mshopts = meshgen('ef',fh,'bou',gdat,'nscreen',5,'plot_on',1,'itmax',50);  
                                                
% now build the mesh with your options and the edge function.
mshopts = mshopts.build; 

%% Get out the msh class from meshgen
m = mshopts.grd;

%plot(m,'tri'); 

%plot(m,'reso');

%% Interpolate on the bathy and gradients (automatically loops over all data)
% m = interp(m,gdat); 
% % ensure max depth in domain is 1 m (we find this step useful for coastal
% % meshes to help with connectivity through narrow channels)
% m.b = max(m.b,1); 
% 
% plot(m,'b');

% Save as a msh class
save('PRVI_msh.mat','m');

% Write an ADCIRC fort.14 compliant file to disk.
%write(m,'PRVI_mesh')

