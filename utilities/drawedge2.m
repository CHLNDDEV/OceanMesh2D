function drawedge2(pp,e2,color)
%DRAWEDGE2 draw EDGE2 elements defined by [PP,E2]. Here, PP
%is an array of positions & E2 is an array of edge-indexing.
if nargin < 3
    color = 'k';
end
%-- draw lines in R^2
xx = [pp(e2(:,1),1), ...
    pp(e2(:,2),1), ...
    NaN * ones(size(e2,1),1)]';
yy = [pp(e2(:,1),2), ...
    pp(e2(:,2),2), ...
    NaN * ones(size(e2,1),1)]';

line('xdata',xx(:), ...
    'ydata',yy(:), ...
    'color',color, ...
    'linewidth',2.0);


end



