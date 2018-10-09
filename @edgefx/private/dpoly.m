function d = dpoly(obj,feat,varargin)
% d = dpoly(obj,feat,bbox)
% obj contains the mesh points and the outer, mainland, inner polygons,
% bbox (the bounding box)
% feat contains inpoly_flip to check whether to flip the inpoly result or not

% d is the distance from point, p to closest point on polygon
% (d is  negative if inside the bounded polygon, pv and positive if outside)
% by Keith Roberts and William Pringle 2017-2018.

%% Doing the distance calc
if ~isempty(obj.lmsl)
    pv = [obj.lmsl.mainland; obj.lmsl.inner];
else
    pv = [feat.mainland; feat.inner];
end
if nargin == 3
    p  = varargin{1};
else
    [xg,yg]=CreateStructGrid(obj);
    p = [xg(:),yg(:)];
    clearvars xg yg
end
pv1 = pv; % dup for inpoly to work
pv1(isnan(pv(:,1)),:) = [];

[~,d] = WrapperForKsearch(pv1, p,1);
% % Do some projection
% lon_mi = min(min(pv1(:,1)),min(p(:,1)));
% lon_ma = max(max(pv1(:,1)),max(p(:,1)));
% lat_mi = min(min(pv1(:,2)),min(p(:,2)));
% lat_ma = max(max(pv1(:,2)),max(p(:,2)));
% m_proj('Transv','long',[lon_mi lon_ma],'lat',[lat_mi lat_ma]);
% [pv1(:,1),pv1(:,2)] = m_ll2xy(pv1(:,1),pv1(:,2));
% p1 = p;
% [p1(:,1),p1(:,2)] = m_ll2xy(p(:,1),p(:,2));
% 
% % Find the nearest neighbors
% idx = WrapperForKsearch(pv1', p1',1);
% % The vector of long lat pairs
% long = zeros(length(p)*2,1);
% lat  = zeros(length(p)*2,1);
% long(1:2:end) = p1(:,1); long(2:2:end) = pv1(idx,1);
% lat(1:2:end) = p1(:,2); lat(2:2:end) = pv1(idx,2);
% % Get spherical earth distances
% d = m_xydist(long,lat); d = d(1:2:end)*1e3;

%% Doing the inpoly check
% If inpoly m file we need to get the edges to pass to it to avoid the
% issues of NaNs
if ~isempty(obj.lmsl)
   edges = Get_poly_edges( [obj.lmsl.outer; obj.lmsl.inner] );
   in = inpoly(p,[obj.lmsl.outer; obj.lmsl.inner],edges);
else
    edges = Get_poly_edges( [feat.outer; feat.inner] );
    [in] = inpoly(p,[feat.outer; feat.inner],edges);
end
% d is negative if inside polygon and vice versa.
if feat.inpoly_flip
    d = (-1).^(~in).*d;
else
    d = (-1).^( in).*d;
end

% % IF OUTSIDE BUT APPEARS INSIDE
% bad = find((p(:,1) < obj.bbox(1,1) | p(:,1) > obj.bbox(1,2) | ...
%             p(:,2) < obj.bbox(2,1) | p(:,2) > obj.bbox(2,2)) & d < 0);
% if nnz(bad) > 0
%     d(bad) = -d(bad);
% end
%EOF
end



