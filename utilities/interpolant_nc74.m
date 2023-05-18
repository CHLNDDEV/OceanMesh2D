function [x_grid_indom,y_grid_indom,wind_grid_indom] = interpolant_nc74(x,y,wind,n_grid)
% function to interpete the original unstrural wind into structural grid
%
%  ref: https://www.mathworks.com/help/matlab/ref/boundary.html
%  ref: https://www.mathworks.com/help/matlab/ref/inpolygon.html
%
%  Author:      Jiangchao Qiu, (MIT/ESSG; email:qiujch24@mit.edu)
%  Created:     May 14 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
[x_grid,y_grid]   = meshgrid(linspace(min(x),max(x),n_grid),linspace(min(y),max(y),n_grid));
wind_interp       = griddata(x,y,wind,x_grid,y_grid);
x_grid_reshape    = reshape(x_grid,[],1);
y_grid_reshape    = reshape(y_grid,[],1);
wind_intp_reshape = reshape(wind_interp,[],1);

%% drop values that are not in the mesh area
k               = boundary(x,y);
x_boundary      = x(k);
y_boundary      = y(k);
in              = inpolygon(x_grid_reshape,y_grid_reshape,x_boundary,y_boundary);
x_grid_indom    = x_grid_reshape(in);
y_grid_indom    = y_grid_reshape(in);
wind_grid_indom = wind_intp_reshape(in);

end