clearvars ; close all; clc ;

minh = 0.01;
maxh = 0.20;
grade = 0.15;

poly =  [0.2483    0.4942
    0.1838    0.6080
    0.2851    0.7861
    0.4741    0.8445
    0.5616    0.8883
    0.6953    0.7277
    0.6745    0.3891
    0.8474    0.3161
    0.9211    0.6168
    0.9764    0.6752
    0.9741    0.1876
    0.7967    0.0241
    0.6838    0.0358
    0.3658    0.1905
    0.3520    0.345
    0.2483    0.4942];

fh = @(x) min(grade*abs(x(:,1)) + minh,maxh);

ID_pfix = []; 
[p,t] = mesh1d(poly,fh,minh,[]);

figure; drawedge2(p,t);
hold on; plot(p(:,1),p(:,2),'s','MarkerFaceColor','r');
axis equal; axis off; drawnow;
