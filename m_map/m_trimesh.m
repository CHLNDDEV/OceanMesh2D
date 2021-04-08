function h = m_trimesh(tri,long,lat,z)
%  M_CONTOURF Adds filled contours to a map
%    M_CONTOUR(LONG,LAT,DATA,...) 
%
global MAP_PROJECTION 

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

[X,Y]=m_ll2xy(long,lat,'clip','on');  

hold on; trimesh(tri,X,Y,z,'facecolor', 'none'); 
%trimesh(tri,X,Y,z,'facecolor', 'flat', 'edgecolor', 'none'); 

%m_coast('patch','red') ; 

end
