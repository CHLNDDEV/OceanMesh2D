function [pts, vals] = edgefx2pts(fh)
% takes an edgefunction from OceanMesh2D and outputs a series of points and
% values for use in fastscatter
% Coleman Blakely
% January 22, 2020
x = fh.F.GridVectors{1};
y = fh.F.GridVectors{2}';
val = fh.F.Values;
vals = reshape(val,[],1);
[lat, lon] = meshgrid(x,y);
c = cat(2,lat',lon');
pts = reshape(c,[],2);
end