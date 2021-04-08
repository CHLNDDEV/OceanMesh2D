function [h]=m_triplot(long,lat,tri)
%  M_CONTOURF Adds filled contours to a map
%    M_CONTOUR(LONG,LAT,DATA,...) 
%
global MAP_PROJECTION 

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end
bcol=[.8,.9,1];
[X,Y]=m_ll2xy(long,lat,'clip','on');  

hold on; trimesh(tri,X,Y,0*X,'facecolor',bcol,'edgecolor','k');
view(2)
%triplot(tri,X,Y,'k','facecolor',bcol); 

%m_coast('patch','red') ; 

end
