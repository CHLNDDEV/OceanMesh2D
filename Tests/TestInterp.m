%% Test for the msh.interp() method using topo-bathy
%%  DEMs with following grid spacings: 
% i)   constant dx = constant dy 
% ii)  constant dx ~= constant dy 
% iii) varying dx ~= varying dy

clearvars; clc;

addpath('..')
addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../m_map/'))

% add the geospatial data filenames
coastline = 'GSHHS_f_L1';
DEMS = ["CUDEM_equal.nc" "CUDEM_unequal.nc" "CUDEM_varying.nc"];

% meshing parameters
min_el    = 40;          % minimum resolution in meters.
max_el    = 500;         % maximum resolution in meters. 
grade     = 0.25;        % mesh grade in decimal percent. 
dis       = grade; 	 % increasing resolution with distance from shoreline.

volumes = [];
% loop over the different DEMs
for dem = DEMS
   dem_name = dem{1}; 
   %x = ncread(dem_name,'lon');
   %y = ncread(dem_name,'lat');
   %min(diff(x))
   %max(diff(x))
   %min(diff(y))
   %max(diff(y))
   %continue  
 
   gdat = geodata('shp',coastline,'dem',dem_name,'h0',min_el);
 
   fh = edgefx('geodata',gdat,...
           'dis',dis,'max_el',max_el,'g',grade);
          
   mshopts = meshgen('ef',fh,'bou',gdat,...
                  'plot_on',0,'proj','utm');
   mshopts = mshopts.build; 

   m = mshopts.grd; 

   m = interp(m,dem_name); 

   %plot(m,'type','bmesh')
   %print([dem_name(1:end-3) '_bmesh.png'],'-r300','-dpng')

   % compute the volume of the mesh
   x = m.p(:,1); y = m.p(:,2);
   areas = polyarea(x(m.t)',y(m.t)');
   [~,bc] = baryc(m);
   volumes = [volumes areas*bc];
end

% Check the different mesh volumes
disp('Mesh volumes:')
disp(volumes)

volume_diff = (max(volumes)-min(volumes))/mean(volumes);
disp('Maximum relative volume difference:')
disp(volume_diff)

if volume_diff > 0.005
    error('Mesh volumes are too disparate with different dems')
    disp('Not Passed: Interp');
else
    disp('Mesh volumes are within 0.5% of each other')
    disp('Passed: Interp');
end
