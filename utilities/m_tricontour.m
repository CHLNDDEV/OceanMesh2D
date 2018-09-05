function [C,h] = m_tricontour(tri,p,z,N)
%  M_CONTOURF Adds filled contours to a map
%    M_CONTOUR(LONG,LAT,DATA,...)
%
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end

[X,Y] = m_ll2xy(p(:,1),p(:,2),'clip','on');
hold on; [C,h] = tricontour([X,Y],tri,z,N);
m_grid('linest','-');
end
