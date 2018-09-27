function [ hh_m ] = ConvertToWGS84( yg, hh_m )
%CONVERTTOWGS84 Convert a grid of points describing edgelength in space in
% planar metres to WGS84 degrees.
% Given a structured grid with edgelength's defined in hh_m,
% convert from planar metres to WGS84 degrees using inverse Haversine
% 
% INPUTS:
% YG : A grid of points in WGS*4 containing the y-locations of these
%      points.
% HH_M : A grid the same size as XG and YG containing the edgelength in
%      planar metres 
% OUTPUTS: 
% HH_M : A grid the same size as XG and YG containing the edglength in
%        WGS84 degrees assuming the edgelength is orientated in the x-direction at
%        each (XG,YG) point. 

% Ensure to remove problem at pole
ul = 89;
yg(yg > ul) = ul; yg(yg < -ul) = -ul;

% Convert to projected coordinates 
Re = 6378.137e3; % <-radius of Earth

% We use a simple inverse Harvesine formula by assuming that the hh_m 
% is applied along a latitude parallel, thus the latitude is constant
% between two points and we only need to solve for the difference in
% longitude (the new hh_m)
el = sin(hh_m./(2*Re));
hh_m = 2*asind(el./cosd(yg));

end

