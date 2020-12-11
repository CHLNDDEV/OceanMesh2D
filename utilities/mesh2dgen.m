function mfp = mesh2dgen( polygon, fh )
% mfp = mesh2dgen( polygon, fh )
% Generates a mesh using mesh2d based on a nan-delimited polygon on wgs84
% lat-lon coordinates and the oceanmesh2D edgefx class, fh
% We make the assumption that floodplain domain is relatively small so that
% projection is not that important..
%
% This section of the code used `mesh2d` by DR. Darren Engwirda
% https://github.com/dengwirda/mesh2d

opts.iter = 100;
opts.kind = 'delaunay';
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
