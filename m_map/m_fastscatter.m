function [h]=m_fastscatter(long,lat,z)
% m_fastscatter uses fastscatter to plot
%
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end

[X,Y]=m_ll2xy(long,lat,'clip','on');
hold on; h=fastscatter(X,Y,z); 
end
