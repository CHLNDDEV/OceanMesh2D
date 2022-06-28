function mfp = mesh2dgen( polygon, fh, kind, iter)
% mfp = mesh2dgen( polygon, fh, kind, iter)
% Generates a mesh using mesh2d based on a nan-delimited polygon on wgs84
% lat-lon coordinates and the oceanmesh2D edgefx class, fh
% We make the assumption that floodplain domain is relatively small so that
% projection is not that important..
% kind: 'delaunay' (default) or 'delfront' method for mesh generation
% iter: maximum allowable iterations (default is 100) 
%
% This section of the code used `mesh2d` by DR. Darren Engwirda
% https://github.com/dengwirda/mesh2d

if nargin < 3 || isempty(kind)
    opts.kind = 'delaunay';
else
    opts.kind = kind;
end
if nargin < 4 || isempty(iter)
    opts.iter = 100;
else
    opts.iter = iter;
end
opts.ref1 = 'preserve';

[node,edge] = getnan2(polygon);

if isa(fh,'edgefx')
    fh.F.Values = fh.F.Values/111e3;
    hfun = @(p)fh.F(p);
    opts.siz1 = 1 + 0.03;
    opts.siz2 = opts.siz1;
else
    hfun = fh;
end

% make the delaunay refinement mesh
[p2,etri,t2,tnum] = refine2_om(node,edge,[],opts,hfun);
% smooth2
[p2,~,t2,~] = smooth2(p2,etri,t2,tnum);

mfp = msh();
mfp.p = p2;
mfp.t = t2;

end
