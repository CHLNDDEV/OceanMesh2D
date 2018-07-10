function [IX,IX1,IX2] = FindLinearIdx(x,y,lon,lat)
% given points in vectors x,y (np x 1) find their linear indices IX (np x 1)
% from a matrix of x-locations X and y-locations Y of both size nx columns and ny rows.
% lon and lat must be matrices created by ndgrid. 
% kjr, 20171210 chl, und.
ny = size(lon,1);
nx = size(lon,2);
np = numel(x);

X = reshape(lon,[],1);
Y = reshape(lat,[],1);

x = x(:);
y = y(:);

dx  = X(2)-X(1);
dy  = dx;

IX1 = (x-X(1))./dx + 1;
IX2 = (y-Y(1))./dy + 1;

IX1 = round(IX1);
IX2 = round(IX2);

IX1 = max([IX1,ones(np,1)],[],2);
IX1 = min([IX1,ny*ones(np,1)],[],2);
%
IX2 = max([IX2,ones(np,1)],[],2);
IX2 = min([IX2,nx*ones(np,1)],[],2);

IX = sub2ind([ny,nx],IX1,IX2);
end