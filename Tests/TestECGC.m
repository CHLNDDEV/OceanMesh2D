% Run example 3 to test multi-scale meshing, bound_courant_number, and interp
clearvars; clc

addpath('..')
addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../m_map/'))

ERR_TOL = 0.05;
ERR2_TOL = 0.01;
QUAL_TOL = 0.25;
Re2 = 111^2;

bbox = [-71.6 42.7; -64 30; -80 24; -85 38; -71.6 42.7]; %polygon boubox
min_el    = 1e3;  		        % minimum resolution in meters.
max_el    = 50e3; 		        % maximum resolution in meters. 
wl        = 30;                 % 60 elements resolve M2 wavelength.
dt        = [50,0.5,6.0];      % Automatically set timestep based on nearshore res
grade     = 0.15;               % mesh grade in decimal percent. 
R         = 3; 			        % Number of elements to resolve feature.
  

dem       = 'SRTM15+.nc';
coastline = 'GSHHS_f_L1';
gdat1 = geodata('shp',coastline,'dem',dem,'h0',min_el,...
                'bbox',bbox);
            
fh1 = edgefx('geodata',gdat1,...
             'fs',R,'wl',wl,'max_el',max_el,...
             'dt',dt,'g',grade);
          
mshopts = meshgen('ef',fh1,'bou',gdat1,...
                  'plot_on',1,'proj','lam');
mshopts = mshopts.build; 

m = mshopts.grd; 

m = interp(m,gdat1); 

%  Set desired courant criteria and time step parameter
dt     = 50; %[s]
max_cr = 4.0; 
min_cr = 0.25;
max_it = 20;
m = bound_courant_number(m,dt,max_cr,min_cr,max_it);

Crnew = CalcCFL(m,dt);

if sum(Crnew > max_cr) > 0
    error(['Unable to bound maximum Courant number. Got ',...
        num2str(max(Crnew)),' expecting < ' num2str(max_cr)]);
else
    disp(['Bounded maximum Courant number. Got ',...
        num2str(max(Crnew)),' expecting < ' num2str(max_cr)]);
end
if sum(Crnew < min_cr) > 0
    error(['Unable to bound minimum Courant number. Got ',...
        num2str(min(Crnew)), ' expecting > ' num2str(min_cr)]);
else
    disp(['Bounded minimum Courant number. Got ',...
        num2str(min(Crnew)), ' expecting > ' num2str(min_cr)]);
end

disp('Passed: ECGC');
