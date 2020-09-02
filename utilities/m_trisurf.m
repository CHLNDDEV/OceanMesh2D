function [h]=m_trisurf(tri,long,lat,z,facecolor)
%  h = m_trisurf(tri,long,lat,z,varargin)
%  
global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
    disp('No Map Projection initialized - call M_PROJ first!');
    return;
end

% facecolor
if nargin < 5
    facecolor = 'interp';
end

% Plotting with a trick to fill in the whole plotting space.
[Xp,Yp] = m_ll2xy(long,lat,'clip','patch');
[X,Y] = m_ll2xy(long,lat,'clip','point');
tx = X(tri); ty = Y(tri);
I = (isnan(tx(:,1)) | isnan(ty(:,1))) & ...
    (isnan(tx(:,2)) | isnan(ty(:,2))) & ...
    (isnan(tx(:,3)) | isnan(ty(:,3)));
hold on; h = trisurf(tri(~I,:),Xp,Yp,z,'facecolor', facecolor, 'edgecolor', 'none');
end