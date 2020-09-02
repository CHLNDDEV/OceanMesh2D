% Run example 5b to build a mesh around Jbay with the weirs
run('../Examples/Example_5b_JBAY_w_weirs.m')

ERR_TOL = 0.05;
ERR2_TOL = 0.01;
QUAL_TOL = 0.25;
Re2 = 111^2;
% p: [23957×2 double]
% t: [42853×3 double]
% Area = 211.19 km^2
% Volume = 2.0782 km^3
if abs(length(m.p) - 23957)/23957 > ERR_TOL
    error(['Incorrect number of points for Example_5b_JBAY_w_weirs.m. Got ',...
        num2str(length(m.p)),' expecting 23957 +- 5% points']);
    exit(1)
end
if abs(length(m.t) - 42853)/42853 > ERR_TOL
    error(['Incorrect number of triangles for Example_5b_JBAY_w_weirs.m. Got ',...
        num2str(length(m.t)),' expecting 42853 +- 5% elements']);
    exit(1)
end
X = m.p(:,1); Y = m.p(:,2);
Area = sum(polyarea(X(m.t)',Y(m.t)').*cosd(mean(Y(m.t)')))*Re2;
if abs(Area - 211.19)/211.19 > ERR2_TOL
    error(['Incorrect area for Example_5b_JBAY_w_weirs.m. Got ',...
        num2str(Area),' km^2 expecting 211.19 km^2 +- 1%']);
    exit(1)
end
bc = (m.b(m.t(:,1))+m.b(m.t(:,2))+m.b(m.t(:,3)))/3;
bc = bc/1e3; % convert to km
Volume = sum(polyarea(X(m.t)',Y(m.t)').*cosd(mean(Y(m.t)')).*bc')*Re2;
if abs(Volume - 2.0782)/2.0782 > ERR2_TOL
    error(['Incorrect volume for Example_5b_JBAY_w_weirs.m. Got ',...
        num2str(Volume),' km^3 expecting 2.0782 km^3 +- 1%']);
    exit(1)
end
if mshopts.qual(end,3) < QUAL_TOL
    error(['Incorrect min. element quality for Example_5b_JBAY_w_weirs.m. Got ',...
        num2str(mshopts.qual(end,3)),' expecting > 0.25']);
    exit(1)
end
if m.bd.nbou ~=2 || isempty(m.bd.ibconn)
    error('Weirs were not passed in Example_5b_JBAY_w_weirs.m correctly.');
    exit(1)
end
disp('Passed: Example_5b_JBAY_w_weirs.m');
