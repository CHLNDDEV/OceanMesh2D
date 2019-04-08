%%
% Example_MeshFloodplain: Mesh a floodplain using five Post-Hurricane Sandy 
% LiDAR 1/9 arc second tiles with shoreline 
% edge locking in a two-step meshing approach. 
% This two-step approach carefully preserve the mesh's 
% shoreline representation in a mesh that extends to a higher geometric countour.
%%%% THIS EXAMPLE TAKES ROUGHLY 25-MINUTES TO COMPLETE %%%%
%%
clearvars; close all; clc;

addpath(genpath('utilities/'))
addpath(genpath('datasets/'))
addpath(genpath('m_map/'))
%% STEP 1: mesh underwater portion of domain first!

% Shoreline dataset, this polygon extends until mean high water outside of
% the local tides' range. 
SHP  = 'us_medium_shoreline_polygon';
DEM  = 'topo15_compressed.nc';

% Construct an outer domain that encapsulates all local tiles. 
BBOX  = [-76 -70; 35 42] ; 
DIS   = 0.15 ; 
GRADE = 0.15 ; 
MIN_EL= 1E3 ; 
MAX_EL= 10E3 ;
SLP   = 20 ; 
WL    = 20 ; 

% Build boundary of outer mesh. 
gdat{1} = geodata('shp',SHP,...
    'h0',MIN_EL,'bbox',BBOX); 

% Build a simple edgefunction for this coarse outer domain.
fh{1}   = edgefx('dis',GRADE,'geodata',gdat{1},...
                 'g',GRADE,'max_el',MAX_EL) ; 

%% Build local high resolution (25-m minimum size) insets 
%% in only areas where LiDAR data is available. 
MIN_EL    = 25 ; % MINIMUM ELEMENT SIZE
MAX_EL_NS = 125 ;% MAXIMUM ELEMENT SIZE NEARSHORE 
MAX_EL    = 1e3; % MAXIMUM ELEMENT SIZE IN THE ENTIRE DOMAIN
FS        = -5 ; % FEATURE SIZE THAT SCALES RESOLUTION WITH SHORELINE WIDTH
GRADE     = 0.15 ; % INTER-ELEMENTAL EXPANSION RATE

% Five 1/9 arc second NCEI tiles from the Post-Sandy dataset referenced to
% MHW. Available from: 
% https://www.ngdc.noaa.gov/mgg/inundation/sandy/sandy_geoc.html
DEMFILES = {'ncei19_n39x25_w075x25_2014v1.nc',...
    'ncei19_n39x25_w075x00_2014v1.nc',...
    'ncei19_n39x25_w074x75_2014v1.nc',...
    'ncei19_n39x00_w075x50_2014v1.nc',...
    'ncei19_n39x00_w075x25_2014v1.nc'};

% Loop over all the tiles of the LiDAR dataset.
for i = 1 : 5
    
gdat{i+1} = geodata('shp', SHP,...
    'dem',DEMFILES{i},...
    'h0',MIN_EL);

fh{i+1} = edgefx('geodata',gdat{i+1},...
    'fs',FS,...
    'max_el_ns',MAX_EL_NS,...
    'max_el',MAX_EL,...
    'g',GRADE);

end
%% CONSTRUCT MESH 
% Note, we set a low disjoint cutoff (dj_cutoff) to retain as much
% hydrualic connectivity of the mesh as possible. A larger value of 
% dj_cutoff will more aggressively simplify the shoreline representation.
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'dj_cutoff',0.001);

mshopts = mshopts.build;

m = mshopts.grd;

plot(m,'tri',0)

% SAVE THE OCEANSIDE MESH AS MSH OBJ TO DISK FOR LATER USE
muw=m;

save OceanSide muw

%% STEP 2: mesh overland while locking shoreline edges and points
%% created from step 1. 

% Only difference in edgefunction is a maximum element size overland of
% 250-m. Note the second line of MAX_EL, with the trailing 
% options maximum and minimum depth above the geoid to enforce 
% the maximum element size constraint.
MAX_EL    = [1e3 0 -inf
             250  inf 0]; 
         
% Loop over the LiDAR tiles once again     
for i = 1:5
    
    % Ensure the feature-size calculation is based on the shoreline
    % datset. 
    gdat_lmsl{i} = geodata('shp',SHP,...
        'dem',DEMFILES{i},...
        'h0',MIN_EL);
    
    % Here we choose to mesh the entire tile. Later on after the mesh is
    % finished and bathymetry is interpolated, we can prune the overland
    % extents according to depth. 
    gdat{i+1} =  geodata('outer', gdat_lmsl{i}.boubox,...
        'dem',DEMFILES{i},...
        'h0',MIN_EL);
    
    gdat{i+1}.inner = [] ;
    
    fh{i+1} = edgefx('geodata',gdat{i+1},...
        'fs',FS,...
        'max_el_ns',MAX_EL_NS,...
        'max_el',MAX_EL,...        
        'lmsl',gdat_lmsl{i},...
        'max_el',MAX_EL,...
        'g',GRADE);
end

%% EXTRACT SHORELINE CONSTRAINTS FROM OCEANSIDE MESH
[pfix,egfix]=muw.extractFixedConstraints; 

% We are going to constrain edges and points only in the high resolution
% insets, which are cell-array entries 2 through 6 in gdat{} and fh{}
% arrays (since box 1 is the outer domain and box 2-6 are the hi-res
% insets).
fixboxes(1) = 0 ;   % don't constrain
fixboxes(2:6) = 1 ; % constrhain

% Pass constraints and fixboxes flags to meshgen.
mshopts = meshgen('ef',fh,'bou',gdat,...
    'pfix',pfix,'egfix',egfix,'fixboxes',fixboxes,...
    'plot_on',1,'nscreen',1);

mshopts = mshopts.build;

mol = mshopts.grd;

%% Quality control: visually inspect mesh!!!!
bou = mol.getBoundaryOfMesh; 

plot(mol,'tri',0); 

hold on; plot(bou(:,1),bou(:,2),'r-','linewi',2); 

%% INTERPOLATE BATHY FROM LiDAR USING SPECIAL PROCEDURE
% This procedure interpolates data underwater FIRST using a gridscale
% averaging, then overland using a larger gridscale averaging to
% produce smoother overland topography. It also ensures the shoreline component is
% at leasts 1-m below sea level. 
mol = interpFP(mol,gdat,muw) ;

plot(mol,'b',0) ;

demcmap([-5 5]);

%%
% At this stage it would appropriate to remove portions of the mesh far
% overlnd out of the reach of a potential flood.
% Here we remove all elements with an average depth greater than 10-m above
% the geoid of the DEMs used. 
mol2 = mol.pruneOverlandMesh(10) ; 

plot(mol2,'tri',0) ; 

plot(mol2,'b',0) ;

demcmap([-5 5]);

% Ensure the CFL is sufficiently small for numerical stability. 
max(CalcCFL(mol2,1))
