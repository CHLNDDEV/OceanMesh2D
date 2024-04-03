function setup_oceanmesh2d()
% add paths for OceanMesh2D

BasePath = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(BasePath,'datasets')))
addpath(genpath(fullfile(BasePath,'utilities')))

% add m_map
if exist('m_map','dir') == 7
    addpath(fullfile(BasePath,'m_map'))
end

addpath(BasePath)
