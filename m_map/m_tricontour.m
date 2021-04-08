function H = m_tricontour(tri,p,z,N,C)
%   H = m_tricontour(tri,p,z,N,C)
%   N-> number of contour levels
%   C-> linespec color
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end

[X,Y] = m_ll2xy(p(:,1),p(:,2),'clip','on');
hold on; 
if nargin < 5
    H = tricontour(tri,X,Y,z,N);
else
    H = tricontour(tri,X,Y,z,N,C);
end

end
