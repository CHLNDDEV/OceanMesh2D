function [ hh_m ] = ConvertToPlanarMetres( xg, yg, hh_m ) 
%CONVERTTOPLANARMETRES
% Given a structured grid with edgelength's defined in hh_m,
% convert it to planar metres using the Haversine formula assuming
% a spherical earth of radius 6378.137km.
% 
% INPUTS:
% XG : A grid of points in WGS84 containing the x-locations of these
%      points.
% YG : A grid of points in WGS*4 containing the y-locations of these
%      points.
% HH_M : A grid the same size as XG and YG containing the edgelength in
%      WGS84 degrees 
% OUTPUTS: 
% HH_M : A grid the same size as XG and YG containing the edglength in
%        planar metres assuming the edgelength is orientated in the x-direction at
%        each (XG,YG) point. 

% Try to remove problem near pole
ul = 89;
yg(yg > ul) = ul; yg(yg < -ul) = -ul;

ptsx = xg(:)'; 
ptsy = yg(:)'; 

ptsx2 = xg(:)'+hh_m(:)'; 
ptsy2 = yg(:)'; 

blahx = [ptsx; ptsx2]; 
blahy = [ptsy; ptsy2]; 

[temp]=m_lldist(blahx(:),blahy(:),2);

hh_m = reshape(temp(1:2:end),size(xg,1),size(xg,2))*1e3; 

end

