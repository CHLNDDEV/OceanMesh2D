function [h]=m_trisurf(tri,long,lat,z)
%  M_CONTOURF Adds filled contours to a map
%    M_CONTOUR(LONG,LAT,DATA,...)
%
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end

[X,Y] = m_ll2xy(long,lat,'clip','on');
hold on; h = trisurf(tri,X,Y,z,'facecolor', 'interp', 'edgecolor', 'none');
end
