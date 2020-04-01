% Run example 1 to build a mesh around NZ testing various basic functionality 
% is working!
run('../Examples/Example_1_NZ.m')

NP_TOL = 500;
NT_TOL = 1500;
QUAL_TOL = 0.15;
% p: [5968×2 double]
% t: [9530×3 double]
% Element qual. metrics
% MEAN      LOWER 3RD  MIN.
% 0.9325    0.7452    0.3548
if abs(length(m.p) - 5968) > NP_TOL
    error(['Incorrect number of points for Example1_NZ.m. Got ',...
        num2str(length(m.p)),' expecting 5968 +- 500 points']);
    exit(1)
end
if abs(length(m.t) - 9530) > NT_TOL
    error(['Incorrect number of triangles for Example1_NZ.m. Got ',...
        num2str(length(m.t)),' expecting 9530 +- 1500 elements']);
    exit(1)
end
if abs(mshopts.qual(end,3) - 0.3548) > QUAL_TOL
    error(['Incorrect min. element quality for Example1_NZ.m. Got ',...
        num2str(mshopts.qual(end,3)),' expecting 0.3548 +- 0.05']);
    exit(1)
end
disp('Passed: Example1_NZ.m');

