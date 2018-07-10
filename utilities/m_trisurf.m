function [h]=m_trisurf(tri,long,lat,z,cmap)
%  M_CONTOURF Adds filled contours to a map
%    M_CONTOUR(LONG,LAT,DATA,...)
%
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end;
m_grid('linest','-');

[X,Y]=m_ll2xy(long,lat,'clip','on');
colormap(cmap); 
hold on; h=trisurf(tri,X,Y,z,'facecolor', 'flat', 'edgecolor', 'none');
shading interp;
end
