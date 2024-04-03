% Run example 1 to build a mesh around NZ testing various basic functionality
% is working!
run('../Examples/Example_1_NZ.m')

PREFIX = 'Example_1_NZ';

NP_TOL = 500;
NT_TOL = 1500;
QUAL_TOL = 0.25; %[mqa for medium cleaning option]
% p: [5968x2 double]
% t: [9530x3 double]
% Element qual. metrics
% MEAN      LOWER 3RD  MIN.
% 0.9325    0.7452    0.3548
TARGET = 5968;
VALUE = length(m.p);
if abs(VALUE - TARGET) > NP_TOL
    error(['Incorrect number of points for %s. ',...
        'Got %i, expecting %i +- %i.'],...
        PREFIX,VALUE,TARGET,NP_TOL);
end

TARGET = 9530;
VALUE = length(m.t); 
if abs(VALUE - TARGET) > NT_TOL
    error(['Incorrect number of triangles for %s. ',...
        'Got %i, expecting %i +- %i.'],...
        PREFIX,VALUE,TARGET,NT_TOL);
end

VALUE = mshopts.qual(end,3);
if VALUE < QUAL_TOL
    error(['Incorrect min. element quality for %s. ',...
        'Got %4.2f, expecting >%4.2f.'],...
        PREFIX,VALUE,QUAL_TOL);
end

fprintf('Passed: %s\n',PREFIX);
