function d = dpoly(obj,feat,varargin)
% INPUTS: 
% obj contains the mesh points and the outer, mainland, inner polygons,
% bbox (the bounding box)
% feat contains inpoly_flip to check whether to flip the inpoly result or not
%
% OUTPUTS: 
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
    [xg,yg] = CreateStructGrid(obj);
    p = [xg(:),yg(:)];
    clearvars xg yg
end
pv1 = pv; % <-dup for inpoly to work
pv1(isnan(pv(:,1)),:) = []; clear pv;
d = 0*p(:,1);
noblks = ceil(length(p)*2*8*1e-9);
blklen = floor(length(p)/noblks);
ns = 1;
for blks = 1:noblks
    if blks == noblks
        ne = length(p); 
    else
        ne = ns + blklen - 1;
    end
    [~,d(ns:ne)] = WrapperForKsearch(pv1, p(ns:ne,:), 1);
    ns = ne + 1;
end

%% Doing the inpoly check
if ~isempty(obj.lmsl)
    edges = Get_poly_edges( [obj.lmsl.outer; obj.lmsl.inner] );
    in = inpoly(p,[obj.lmsl.outer; obj.lmsl.inner],edges);
else
    edges = Get_poly_edges( [feat.outer; feat.inner] );
    in = inpoly(p,[feat.outer; feat.inner],edges);
end
% d is negative if inside polygon and vice versa.
if feat.inpoly_flip
    d = (-1).^(~in).*d;
else
    d = (-1).^( in).*d;
end
end



