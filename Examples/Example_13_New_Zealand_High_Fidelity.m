% Example_13_NZ_high_fidelity: High-Fidelity Mesh Generation for the South Island of New Zealand
% This script demonstrates the application of multiscale meshing techniques with a focus on high-fidelity options.
% High-fidelity meshing is particularly crucial for areas within a larger regional model that contain engineered structures
% or coastlines. Such areas require precise representation in the mesh connectivity to ensure accurate simulation outcomes.
% This approach enables detailed analysis and modeling of geophysical processes near these critical features, enhancing the
% reliability of simulations in coastal engineering, environmental assessment, and resource management applications.
clearvars; clc;
clearvars; clc;

% Adding directories to the MATLAB path for utilities, datasets, and mapping tools.
addpath('..');
addpath(genpath('../utilities/'));
addpath(genpath('../datasets/'));
addpath(genpath('../m_map/'));

%% STEP 1: Define mesh boundaries and resolution parameters.
% Bounding box coordinates (longitude and latitude) define the area of interest.
bbox1 = [166, 176; -48, -40]; % lon_min, lon_max; lat_min, lat_max
min_el1 = 1e3; % Minimum element size (m).
max_el1 = 100e3; % Maximum element size (m).
max_el_ns1 = 5e3; % Maximum nearshore element size (m).
grade1 = 0.15; % Mesh grading factor.
R1 = 3; % Feature resolution factor.

%% STEP 2: Initialize geographic data for mesh generation.
% Specify coastline dataset for mesh boundary definition.
coastline = 'GSHHS_f_L1';
gdat{1} = geodata('shp', coastline, 'bbox', bbox1, 'h0', min_el1);

%% STEP 3: Generate edge functions for mesh refinement.
% Edge function defines transition of mesh sizes across the domain.
fh{1} = edgefx('geodata', gdat{1}, 'fs', R1, 'max_el_ns', max_el_ns1, 'max_el', max_el1, 'g', grade1);

%% STEP 4: Define a refined mesh area.
% Additional bounding box for targeted refinement (e.g., a specific bay or coastal region).
% In this case Christchurch Bay with a more detailed shoreline vector
% file that contains hardened shoreline features.
coastline2 = 'coastline_epsg4326';


bbox2= [172.61348324,172.85932924;
    -43.67102284,-43.51383230];

min_el2 = 50.0; % Refined minimum element size (m).
max_el_ns2 = 250.0; % Refined maximum nearshore element size (m).
max_el2 = 500; % Refined maximum element size (m).
grade2 = 0.05; % Mesh grading factor for refined area.
R2 = 5; % Feature resolution factor for refined area.

% High-fidelity option for enhanced mesh detail in the refined area.
gdat{2} = geodata('shp', coastline2, 'bbox', bbox2, 'h0', min_el2, 'high_fidelity', 1);
fh{2} = edgefx('geodata', gdat{2}, 'fs', R2, 'max_el_ns', max_el_ns2, 'max_el', max_el2, 'g', grade2);

%% STEP 5: Mesh generation with defined edge functions.
% Configure mesh generation options, including visualization during the process.
mshopts = meshgen('ef', fh, 'bou', gdat, 'plot_on', 1, 'nscreen', 5);

% You can extract the developed constraints
pfix = mshopts.pfix;
egfix = mshopts.egfix;

mshopts = mshopts.build(); % Execute mesh building with specified options.

%% STEP 6: Visualize and export the generated mesh.
% Extract and visualize the generated mesh, highlighting refined areas.
m = mshopts.grd;
m.plot('type', 'tri', 'proj', 'lamb'); % Plot entire mesh.

m.plot('type', 'tri', 'proj', 'none', 'subdomain', bbox2); % Highlight refined area.
hold on; drawedge2(pfix, egfix, 'r'); % Overlay edge constraints.

% Optional: Export mesh to a .2dm file for editing connectivity in GIS software like QGIS.
write(m, 'South_Island_NZ_w_constraints.2dm', '2dm');
