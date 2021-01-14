% Example_12_Riverine_flow: Guidance for creating a riverine flow input file 
% otherwise referred to as a fort.20
%
% You should prepare a column delimited file (csv) that stores the total volume flow
% (Q) (m^3/s) time series of a cross-section where the riverine boundaries in the 
% mesh are located. 
%
% The columns of the csv file should be organized as follows:  
% year,month,day,hour,minute,second,volume_flow_1,volume_flow_2... 
% The column order of the total volume flow (volume_flow_1,volume_flow_2...) 
% must be specified in the order in which the riverine boundaries appear in 
% the fort.14 file, or in which you make the riverine boundaries with the data 
% cursor method in msh.make_bc.
%
% A 'test_make_f20.csv' file has been placed in the Tests folder for reference.
% The volume flow of this file is an hourly average time series, thus DT equals 
% to 3600 (s), if your volume flow is a daily average time series, DT should be  
% set to 86400 (s) rather than 3600 (s).
%
% More details can be found in the two functions of Make_f20_volume_flow.m and 
% Riverflux_distribution.m.
%
% Author:      Jiangchao Qiu                                
% Created:     January 7 2021                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('utilities/'));
addpath(genpath('datasets/'));
addpath(genpath('m_map/'));

%% STEP 1: set mesh extents and set parameters for mesh. 
%% The Pearl River Estuary region
bbox = [112 20; 112 24; 116 24; 116 20; 112 20]; %polygon boubox
min_el    = 100;  		       % minimum resolution in meters.
max_el    = 10e3; 		       % maximum resolution in meters. 
dt        = 5;                 % Automatically set timestep based on nearshore res
grade     = 0.3;               % mesh grade in decimal percent. 
R         = 3; 			       % Number of elements to resolve feature.
  
%% STEP 2: specify geographical datasets and process the geographical data
%% to be used later with other OceanMesh classes...
dem       = 'SRTM15+V2.1.nc';
coastline = 'GSHHS_f_L1';
gdat = geodata('shp',coastline,'dem',dem,'h0',min_el,'bbox',bbox);
%plot(gdat,'shp'); 
%plot(gdat,'dem');  
       
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,'fs',R,'max_el',max_el,'dt',dt,'g',grade);

%% STEP 4: Pass your edgefx class object along with some meshing options and... 
%% build the mesh         
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'proj','lam');
mshopts = mshopts.build; 

%% STEP 5: Plot and save the msh class object
m = mshopts.grd; % get out the msh object
plot(m,'type','tri');

m = interp(m,gdat,'mindepth',1); % interpolate bathy to the mesh with minimum depth of 1 m
m = lim_bathy_slope(m,0.1,0);
CFL = CalcCFL(m,dt); 
max(CFL)
min(CFL)
m = BoundCr(m,dt);

%% make the nodestring boundary conditions
m = make_bc(m,'auto',gdat); 

%% use data cursor to identify vstart and vend of each riverine boundary
m = make_bc(m,'outer',1);% 1,2206,2107,1(flux),22(River)
m = make_bc(m,'outer',1);% 1,14784,14779,1(flux),22(River)
m = make_bc(m,'outer',1);% 1,763,704,1(flux),22(River)

plot(m,'type','bd');  % plot mesh on native projection with boundary conditions
plot(m,'type','b');   % plot bathy on native projection

ts = '01-Jan-2005 00:00'; % start time of the total volume flow time series
te = '01-Jan-2005 23:00'; % end time of the total volume flow time series
DT = 3600; % time step of the total volume flow time series,
           % DT=3600 for hourly data series while DT=86400 for daily data series.
m = Make_f20_volume_flow(m,'test_make_f20.csv',ts,te,DT);

save('Riverine.mat','m'); write(m,'Riverine','20');



