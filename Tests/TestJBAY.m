% Run example 5b to build a mesh around Jbay with the weirs

run('../Examples/Example_5b_JBAY_w_weirs.m')

PREFIX = 'Example_5b_JBAY_w_weirs';

ERR_TOL = 0.05;
ERR2_TOL = 0.01;
QUAL_TOL = 0.25;
Re2 = 111^2;
% p: [25190x2 double]
% t: [45270x3 double]
% Area = 211 km^2
% Volume = 2.075 km^3
TARGET = 25190;
VALUE = length(m.p);
if abs(VALUE - TARGET)/TARGET > ERR_TOL
    error(['Incorrect number of points for %s. ',...
        'Got %i, expecting %i +- %4.2f%%.'],...
        PREFIX,VALUE,TARGET,ERR_TOL*100);
end
TARGET = 45270;
VALUE = length(m.t); 
if abs(VALUE - TARGET)/TARGET > ERR_TOL
    error(['Incorrect number of triangles for %s. ',...
        'Got %i, expecting %i +- %4.2f%%.'],...
        PREFIX,VALUE,TARGET,ERR_TOL*100);
end

TARGET = 211;
X = m.p(:,1); Y = m.p(:,2);
Area = sum(polyarea(X(m.t)',Y(m.t)').*cosd(mean(Y(m.t)')))*Re2;
if abs(Area - TARGET)/TARGET > ERR2_TOL
    error(['Incorrect area for %s. ',...
        'Got %4.2f km^2, expecting %4.2f km^2 +- %4.2f%%.'],...
        PREFIX,Area,TARGET,ERR2_TOL*100);
end

TARGET = 2.075;
bc = (m.b(m.t(:,1))+m.b(m.t(:,2))+m.b(m.t(:,3)))/3;
bc = bc/1e3; % convert to km
Volume = sum(polyarea(X(m.t)',Y(m.t)').*cosd(mean(Y(m.t)')).*bc')*Re2;
if abs(Volume - TARGET)/TARGET > ERR2_TOL
    error(['Incorrect volume for %s. ',...
        'Got %4.2f km^3, expecting %4.2f km^3 +- %4.2f%%.'],...
        PREFIX,Volume,TARGET,ERR2_TOL*100);
end

VALUE = mshopts.qual(end,3);
if VALUE < QUAL_TOL
    error(['Incorrect min. element quality for %s. ',...
        'Got %4.2f, expecting >%4.2f.'],...
        PREFIX,VALUE,QUAL_TOL);
end

if m.bd.nbou ~=2 || isempty(m.bd.ibconn)
    error('Weirs were not passed in %s correctly.',PREFIX);
end

fprintf('Passed: %s\n',PREFIX);
