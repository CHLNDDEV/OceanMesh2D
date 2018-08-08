function [h]=m_fastscatter(long,lat,z,cmap)
% m_fastscatter uses fastscatter to plot
%
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end
m_grid('linest','-');

[X,Y]=m_ll2xy(long,lat,'clip','on');
colormap(cmap); 
hold on; h=fastscatter(X,Y,z); 
end
