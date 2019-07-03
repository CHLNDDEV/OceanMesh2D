function drawedge2(pp,e2)
%DRAWEDGE2 draw EDGE2 elements defined by [PP,E2]. Here, PP
%is an array of positions & E2 is an array of edge-indexing.

ec = [0 0 0 ];

%-- draw lines in R^2
xx = [pp(e2(:,1),1), ...
    pp(e2(:,2),1), ...
    NaN * ones(size(e2,1),1)]';
yy = [pp(e2(:,1),2), ...
    pp(e2(:,2),2), ...
    NaN * ones(size(e2,1),1)]';

line('xdata',xx(:), ...
    'ydata',yy(:), ...
    'color',ec, ...
    'linewidth',2.0);


end



