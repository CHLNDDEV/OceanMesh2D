function [longMAG,latMAG,phiVecMAG,thetaVecMAG]=m_geo2mag(longGEO,latGEO,phiVecGEO,thetaVecGEO)
% M_GEO2MAG  Converts geographic to geomagnetic coordinates.
%   [longMAG,latMAG]=M_GEO2MAG(longGEO,latGEO) converts magnetic 
%   (dipole) coordinates to geographic coordinates. All in units of 
%   degrees with + longitudes east. All variables can be scalar or matrix 
%   but must have the same size. The geomagnetic coordinate system must
%   previously have been initialized with a call to M_COORD.
%
%   Vector rotations can be carried using
% 
%  [lonMAGG,latMAG,phiVecMAG,thetaVecMAG]=M_GEO2MAG(longGEO,latGEO,phiVecGEO,thetaVecGEO)
%
%   where
%
% phiVecGEO   - east component of the vector in geographic coordinates
% thetaVecGEO - north component of the vector in geographic coordinates
% phiVecMAG   - east component of the vector in geomagnetic coordinates
% thetaVecMAG - north component of the vector in geomagnetic coordinates
%
% See also M_COORD, M_MAG2GEO

% References:
%
% Hapgood, M.A., Space Physics Coordinate Transformations:
% A User Guide, Planet. Space Sci., Vol. 40, N0. 5, 1992.

% R. Pawlowicz (rich@ocgy.ubc.ca)
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% Have to have initialized a map first

global MAP_COORDS

if isempty(MAP_COORDS)
  error('No coordinate initialized - call M_PROJ or M_COORD first!');
elseif strcmp(MAP_COORDS.name.name,'geographic')
  error('No Geomagnetic coordinate system initialized - call M_COORD to choose one');
end

if nargin==2
   [longMAG,latMAG]=mc_coords('geo2mag',longGEO,latGEO);
elseif nargin==4
   [longMAG,latMAG,phiVecMAG,thetaVecMAG]=mc_coords('geo2mag',longGEO,latGEO,phiVecGEO,thetaVecGEO);
else
   error('Wrong number of input parameters');
end  

