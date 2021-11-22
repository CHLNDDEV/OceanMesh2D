function [IX,IX1,IX2] = FindLinearIdx(x,y,lon,lat)
% given points in vectors x,y (np x 1) find their linear indices IX (np x 1)
% from a matrix of x-locations X and y-locations Y of both size nx columns and ny rows.
% lon and lat must be matrices created by ndgrid. 
% kjr, 20171210 chl, und.
% wjp, 20201120 making sure dx and dy can be different
%      20200325 implementing method for irregular grid, cleaning up
ny = size(lon,1);
nx = size(lon,2);
np = numel(x);

if nx == 1 && ny == 1
   IX = 1; IX1 = 1; IX2 = 1;
   return
end

% make sure entry points are column vectors
x = x(:);
y = y(:);

% Get the grid spacing in x and y directions 
% (trying both directions so could be meshgrid or ndgrid format)
dx  = diff(lon(:,1));
dy  = diff(lat(1,:));
dx_cutoff = 0.1/111e3; % approx 10 cm 
if (max(dx) - min(dx)) > dx_cutoff || (max(dy) - min(dy)) > dx_cutoff
    % % IRREGULAR GRID (SLOWER)
    
    % convert ndgrid to vector
    lonV = reshape(lon(:,1),1,[]); 
    % repeat the vector along the length of x
    LON = repmat(lonV,np,1);
    % find the closest lon to each point of x
    [~,IX1] = min(abs(x - LON),[],2);
    
    % convert ndgrid to vector
    latV = lat(1,:); 
    % repeat the vector along the length of y
    LAT = repmat(latV,np,1);
    % find the closest lon to each point of y
    [~,IX2] = min(abs(y - LAT),[],2);

else
    % % REGULAR GRID (FASTER)
    dx = dx(1); dy = dy(1);

    IX1 = (x-lon(1,1))/dx + 1;
    IX2 = (y-lat(1,1))/dy + 1;

    IX1 = round(IX1);
    IX2 = round(IX2);

    IX1 = max([IX1,ones(np,1)],[],2);
    IX1 = min([IX1,ny*ones(np,1)],[],2);

    IX2 = max([IX2,ones(np,1)],[],2);
    IX2 = min([IX2,nx*ones(np,1)],[],2);

end

IX = sub2ind([ny,nx],IX1,IX2);

end
