% Test mesh1d (the 1D boundary resampler used internally by meshgen's
% breakline/high-fidelity constraint generation) on a synthetic polygon
% with a linearly-graded sizing function.

clearvars ; close all; clc ;

run('../setup_oceanmesh2d.m')

PREFIX = 'test_1d_original';

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
[p,t,converged] = meshgen.call_mesh1d(poly,fh,minh,[]);

if converged ~= 1
    error('mesh1d did not converge for %s.',PREFIX);
end

NP_TOL = 5;
TARGET = 138;
VALUE = size(p,1);
if abs(VALUE - TARGET) > NP_TOL
    error(['Incorrect number of points for %s. ',...
        'Got %i, expecting %i +- %i.'],...
        PREFIX,VALUE,TARGET,NP_TOL);
end

NT_TOL = 5;
TARGET = 137;
VALUE = size(t,1);
if abs(VALUE - TARGET) > NT_TOL
    error(['Incorrect number of edges for %s. ',...
        'Got %i, expecting %i +- %i.'],...
        PREFIX,VALUE,TARGET,NT_TOL);
end

d = sqrt(sum((p(t(:,1),:)-p(t(:,2),:)).^2,2));

MIN_EDGE_TOL = 0.005; % slack below the target minimum edge length
if min(d) < minh - MIN_EDGE_TOL
    error(['Edge length below minimum resolution for %s. ',...
        'Got %6.4f, expecting >= %6.4f.'],...
        PREFIX,min(d),minh-MIN_EDGE_TOL);
end

MAX_EDGE_TOL = 0.02; % slack above the target maximum edge length
if max(d) > maxh + MAX_EDGE_TOL
    error(['Edge length above maximum resolution for %s. ',...
        'Got %6.4f, expecting <= %6.4f.'],...
        PREFIX,max(d),maxh+MAX_EDGE_TOL);
end

figure; drawedge2(p,t);
hold on; plot(p(:,1),p(:,2),'s','MarkerFaceColor','r');
axis equal; axis off; drawnow;

fprintf('Passed: %s\n',PREFIX);
