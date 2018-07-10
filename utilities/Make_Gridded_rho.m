function [t,lon,lat,p,rho,zeta,dpx,dpy] = Make_Gridded_rho(filename,names)
% function [t,lon,lat,z,rho,zeta,dpx,dpy] = Make_Gridded_rho(filename,names)
%  Make_Gridded_rho                                                 
%  Description: Read NETCDF of grided salinity and temperature. 
%               Use TEOS toolbox to compute profiles of the in-situ density
%  Inputs:      NETCDF file of salinity and temperature from e.g. a 3D
%               ocean model and cell/string array of variable names in 
%               order of time, depths, lon, lat, salinity and tempererature
%  Outputs:     time, lon, lat, depth contours (z) and in-situ density (rho)   
%  Requires:    TEOS toolbox http://www.teos-10.org/software.htm  
%  Author:      William Pringle                                   
%  Created:     Mar 26 2018, Updated May 25 2018 to give back the 
%               baroclinic pressure gradients                        
%                                                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = filename;
t   = ncread(file,names{1});
units = ncreadatt(file,names{1},'units');
C = strsplit(units);
t_o = datetime(strjoin(C(3:4)));
if strcmp(C{1},'hours')
    t = t/24;
elseif strcmp(C{1},'seconds')
    t = t/24/3600;
end
t = t + t_o;
p   = ncread(file,names{2});
lon = ncread(file,names{3});
lat = ncread(file,names{4});
SP   = ncread(file,names{5});
T    = ncread(file,names{6});
zeta   = ncread(file,names{7});
L = length(lat);
N = length(lon);
rho0 = 1000;

% Make lon lat nd-grid vectors and reshape
[lonN,latN] = ndgrid(lon,lat);
SP = reshape(SP,[],length(p));
T = reshape(T,[],length(p));
lonN = reshape(lonN,[],1);
latN = reshape(latN,[],1);
% Get the absolute salinities and conservative temperatures;
SA = gsw_SA_from_SP(SP,p,lonN,latN);
CT = gsw_CT_from_t(SA,T,p);
% Get the density 
rho = gsw_rho(SA,CT,p);
rho = reshape(rho,N,L,[]);

% Let's remove outliers in density


% % Now calculate the pressures
%zeta = reshape(zeta,[],1);
%b = reshape(b,[],1);
% ph = 0*rho;
% % Get the surface pressure anomoly ps = g*rho*zeta
% ph(:,1) = zeta.*rho(:,1);
% 
% % Get the hydrostatic pressure by integrating density
% for kk = 2:length(p)
%     ph(:,kk) = trapz(p(1:kk),rho(:,1:kk),2);
%     % For areas with depth larger than current depth
% end
% ph(:,2:end) = ph(:,1) + ph(:,2:end);
% ph = g*reshape(ph,N,L,[])/rho0;

% The method below of integrating the gradients (instead of integrating to 
% get pressure then taking gradient) is suggested in the NEMO user guide to 
% reduce truncation error
dpx = 0*rho; dpy = 0*rho;

% Get the average dx and dy
dx = m_idist(lon(floor(end/2)),mean(lat),lon(floor(end/2)+1),mean(lat));
dy = m_idist(mean(lon),lat(floor(end/2)),mean(lon),lat(floor(end/2)+1));
bottom = 0*zeta;
for kk = 1:length(p)
    if kk == 1
        % Get the surface pressure anomaly ps = g*rho*zeta
        ps = zeta.*rho(:,:,1);
        %[dpy(:,:,1),dpx(:,:,1)] = my_gradient(ps,dy,dx);
        [dpy(:,:,1),dpx(:,:,1)] = gradient(ps,dy,dx);
%         figure;
%         pcolor(lon,lat,hypot(dpx(:,:,kk),dpy(:,:,kk))')
%         shading interp
%         caxis([0 1e-4])
%         title(p(kk))
    else
%         if kk == 17
%             I = find(lon > -66.35 & lon < -66.25)
%             J = find(lat > 17.73 & lat < 17.77)
%         end
        % Get gradient at each level
        %[dpy(:,:,kk),dpx(:,:,kk)] = my_gradient(rho(:,:,kk)-rho0,dy,dx);
        [dpy(:,:,kk),dpx(:,:,kk)] = gradient(rho(:,:,kk)-rho0,dy,dx);
    end
    % Check where the bottom is
    bottom(~isnan(rho(:,:,kk))) = kk;
end
% Removing areas where there is large gradient in depth
% gb = gradient(bottom,1);
dpy = reshape(dpy,[],length(p));
dpx = reshape(dpx,[],length(p));
% large_gb = 2;
% large_grad = 0.1;
% dptemp = dpy(abs(gb) > large_gb,:);
% dptemp(abs(dptemp)*dx > large_grad) = 0; 
% dpy(abs(gb) > large_gb,:) = dptemp;
% dptemp = dpx(abs(gb) > large_gb,:);
% dptemp(abs(dptemp)*dy > large_grad) = 0; 
% dpx(abs(gb) > large_gb,:) = dptemp;
% Just remove really large gradients
%dpy(abs(dpy)*dy > 0.1) = 0;
%dpx(abs(dpx)*dx > 0.1) = 0;
% Do the integration in backwards order
for kk = length(p):-1:2
    dpx(:,kk) = trapz(p(1:kk),dpx(:,1:kk),2);
    dpy(:,kk) = trapz(p(1:kk),dpy(:,1:kk),2);
    
%     figure;
%     pcolor(lon,lat,hypot(dpx(:,:,kk),dpy(:,:,kk))')
%     shading interp
%     caxis([0 1e-4])
%     title(p(kk))
end
% Reshape back
dpx = reshape(dpx,N,L,[]);
dpy = reshape(dpy,N,L,[]);
%EOF
end

function [rho,pm] = ADCIRC_EqnState2(z,S,T,zeta)
    g = 9.807; 
    Q1 = 196637339; Q2 = 196668.928;
    a1 = 1446045.44; a2 = -10769.1980; a3 = 85.2955498;
    b1 = 579854.265; b2 = -1477.31528; b3 = 419.489109;
    c1 = 2001.22349; c2 = 0.0213851025; c3 = 0.968784141;
    c4 = -0.00654785602; c5 = -0.00000250468726; c6 = 0.0000000628902345;
    c7 = 0.00000000282414763; d4 = 0.0000398630534;
    d1 = 1433.02205; d2 = -9.09231525; d3 = 0.0791654429; 
    e1 = 417.831720; e2 = -1.87581316; e3 = -0.0000387902837;
    f1 = 1.00765828; f2 = 0.000312912597;
    T_min = 0; T_max = 40;
    S_min = 0; S_max = 42;
    S( S < S_min) = S_min; S( S > S_max) = S_max;
    T( T < T_min) = T_min; T( T > T_max) = T_max;
    
    %z = repmat(z,size(S,1),size(S,2),1);
    rho = 0*S; pm = 0*S;
    for kk = 1:length(z)
        if kk == 1
            p = 0*zeta;
            pm(:,:,kk) = p;
            P1 = Q1 + (a1 + (a2 + a3*T(:,:,kk)).*T(:,:,kk)).*...
                 T(:,:,kk) + (b1 + b2*T(:,:,kk) + b3*S(:,:,kk)).*S(:,:,kk);
            P2 = Q2 + (d1 + (d2 + (d3 + d4*T(:,:,kk)).*T(:,:,kk)).*...
                       T(:,:,kk)).*T(:,:,kk) + ...
                      (e1 + (f1 + f2*T(:,:,kk).*T(:,:,kk)).* ...
                       sqrt(S(:,:,kk)) + (e2 + e3*T(:,:,kk).* ... 
                       T(:,:,kk)).*T(:,:,kk)).*S(:,:,kk);
            rho(:,:,kk) = P1./P2;
        else
            p = p + (z(kk)-z(kk-1))*g*1d-4*rho(:,:,kk-1);
            pm(:,:,kk) = p;
            P1 = Q1 + (a1 + (a2 + a3*T(:,:,kk)).*T(:,:,kk)).* ...
                 T(:,:,kk) + (b1 + b2*T(:,:,kk) + ...
                 b3.*S(:,:,kk)).*S(:,:,kk) + (c1 + c2* ...
                 T(:,:,kk).*T(:,:,kk) + c3*S(:,:,kk) + ...
                 (c4 + ((c5+c6*T(:,:,kk)).*T(:,:,kk) + ...
                 c7*p).*T(:,:,kk)).*p).*p;
            P2 = Q2 + (d1 + (d2 + (d3 + d4*T(:,:,kk)).* ...
                 T(:,:,kk)).*T(:,:,kk)).*T(:,:,kk) + ...
                 (e1 + (f1 + f2*T(:,:,kk).*T(:,:,kk)).* ...
                 sqrt(S(:,:,kk)) + (e2 + e3*T(:,:,kk).* ...
                 T(:,:,kk)).*T(:,:,kk)).*S(:,:,kk);    
            rho(:,:,kk) = P1./(P2 + p);
        end
    end
end
   
function varargout = my_gradient(f,varargin)
%GRADIENT Approximate gradient.
%   [FX,FY] = GRADIENT(F) returns the numerical gradient of the
%   matrix F. FX corresponds to dF/dx, the differences in x (horizontal) 
%   direction. FY corresponds to dF/dy, the differences in y (vertical) 
%   direction. The spacing between points in each direction is assumed to 
%   be one. When F is a vector, DF = GRADIENT(F) is the 1-D gradient.
%
%   [FX,FY] = GRADIENT(F,H), where H is a scalar, uses H as the
%   spacing between points in each direction.
%
%   [FX,FY] = GRADIENT(F,HX,HY), when F is 2-D, uses the spacing
%   specified by HX and HY. HX and HY can either be scalars to specify
%   the spacing between coordinates or vectors to specify the
%   coordinates of the points.  If HX and HY are vectors, their length
%   must match the corresponding dimension of F.
%
%   [FX,FY,FZ] = GRADIENT(F), when F is a 3-D array, returns the
%   numerical gradient of F. FZ corresponds to dF/dz, the differences
%   in the z direction. GRADIENT(F,H), where H is a scalar, 
%   uses H as the spacing between points in each direction.
%
%   [FX,FY,FZ] = GRADIENT(F,HX,HY,HZ) uses the spacing given by
%   HX, HY, HZ. 
%
%   [FX,FY,FZ,...] = GRADIENT(F,...) extends similarly when F is N-D
%   and must be invoked with N outputs and either 2 or N+1 inputs.
%
%   Note: The first output FX is always the gradient along the 2nd
%   dimension of F, going across columns.  The second output FY is always
%   the gradient along the 1st dimension of F, going across rows.  For the
%   third output FZ and the outputs that follow, the Nth output is the
%   gradient along the Nth dimension of F.
%
%   Examples:
%       [x,y] = meshgrid(-2:.2:2, -2:.2:2);
%       z = x .* exp(-x.^2 - y.^2);
%       [px,py] = gradient(z,.2,.2);
%       contour(z), hold on, quiver(px,py), hold off
%
%   Class support for input F:
%      float: double, single
%
%   See also DIFF, DEL2.

%   Copyright 1984-2015 The MathWorks, Inc.

[f,ndim,loc,rflag] = parse_inputs(f,varargin);
nargoutchk(0,ndim);

% bad multiple 
bm = 2;

% Loop over each dimension. 

varargout = cell(1,ndim);
siz = size(f);
% first dimension 
g  = zeros(size(f),class(f)); % case of singleton dimension
h = loc{1}(:); 
n = siz(1);
% Take forward differences on left and right edges
if n > 1
   g(1,:) = (f(2,:) - f(1,:))/(h(2)-h(1));
   g(n,:) = (f(n,:) - f(n-1,:))/(h(end)-h(end-1));
end

% Take centered differences on interior points
if n > 2
   gc = (f(3:n,:)-f(1:n-2,:)) ./ (h(3:n) - h(1:n-2));
end

% Let's check quality
% Get forward and backward differences
gf = (f(3:n,:) - f(2:n-1,:)) ./ (h(3:n)-h(2:n-1));
gb = (f(2:n-1,:) - f(1:n-2,:)) ./ (h(2:n-1)-h(1:n-2));
% Try fill NaNed areas with forward or backward difference
gc(isnan(gc)) = gf(isnan(gc)); gc(isnan(gc)) = gb(isnan(gc)); 
% Fill Really large differences with smaller of forward or backward
bad1 = abs(gf) > bm*abs(gb); 
bad2 =  abs(gb) > bm*abs(gf);
gc(bad1) = gf(bad1);
gc(bad2) = gb(bad2);
% Put centre differences in full gradient
g(2:n-1,:) = gc;

varargout{1} = g;

% second dimensions and beyond
if ndim == 2
    % special case 2-D matrices to support sparse matrices,
    % which lack support for N-D operations including reshape
    % and indexing
    n = siz(2);
    h = reshape(loc{2},1,[]);
    g = zeros(size(f),class(f));
    
    % Take forward differences on left and right edges
    if n > 1
        g(:,1) = (f(:,2) - f(:,1))/(h(2)-h(1));
        g(:,n) = (f(:,n) - f(:,n-1))/(h(end)-h(end-1));
    end
    
    % Take centered differences on interior points
    if n > 2
        gc = (f(:,3:n) - f(:,1:n-2)) ./ (h(3:n) - h(1:n-2));
    end
    
    % Let's check quality
    % Get forward and backward differences
    gf = (f(:,3:n) - f(:,2:n-1)) ./ (h(3:n)-h(2:n-1));
    gb = (f(:,2:n-1) - f(:,1:n-2)) ./ (h(2:n-1)-h(1:n-2));
    % Try fill NaNed areas with forward or backward difference
    gc(isnan(gc)) = gf(isnan(gc)); gc(isnan(gc)) = gb(isnan(gc)); 
    % Fill Really large differences with smaller of forward or backward
    bad1 = abs(gf) > bm*abs(gb); 
    bad2 =  abs(gb) > bm*abs(gf);
    gc(bad1) = gb(bad1);
    gc(bad2) = gf(bad2);
    % Put centre differences in full gradient
    g(:,2:n-1) = gc;
    
    varargout{2} = g; 
end

% Swap 1 and 2 since x is the second dimension and y is the first.
if ndim > 1
    varargout(2:-1:1) = varargout(1:2);
elseif rflag
    varargout{1} = varargout{1}.';
end


%-------------------------------------------------------
function [f,ndim,loc,rflag] = parse_inputs(f,v)
%PARSE_INPUTS
%   [ERR,F,LOC,RFLAG] = PARSE_INPUTS(F,V) returns the spacing
%   LOC along the x,y,z,... directions and a row vector
%   flag RFLAG. 

loc = {};

% Flag vector case and row vector case.
ndim = ndims(f);
vflag = false;
rflag = false;
if isvector(f)
    ndim = 1;
    vflag = true;
    if isrow(f) % Treat row vector as a column vector
        rflag = true;
        f = f.';
    end
end

indx = size(f);

% Default step sizes: hx = hy = hz = 1
if isempty(v)
    % gradient(f)
    loc = cell(1, ndims(f));
    for k = 1:ndims(f)
        loc(k) = {1:indx(k)};
    end
elseif isscalar(v) % gradient(f,h)
    % Expand scalar step size
    if isscalar(v{1})
        loc = cell(1, ndims(f));
        for k = 1:ndims(f)
            h = v{1};
            loc(k) = {h*(1:indx(k))};
        end
        % Check for vector case
    elseif vflag
        loc(1) = v(1);
    else
        error(message('MATLAB:gradient:InvalidInputs'));
    end
elseif ndims(f) == numel(v)  % gradient(f,hx,hy,hz,...)
    % Swap 1 and 2 since x is the second dimension and y is the first.
    loc = v;
    if ndim > 1
        loc(2:-1:1) = loc(1:2);
    end
    % replace any scalar step-size with corresponding position vector
    for k = 1:ndims(f)
        if isscalar(loc{k})
            loc{k} = loc{k}*(1:indx(k));
        end
    end 
else
    error(message('MATLAB:gradient:InvalidInputs'));
end
end
end
% Loop over data in latitudinal direction because GSW are 
% functions of latitude

%SA_sp = cell(L,1);
%CT_sp = cell(L,1);

% [~, pADC] = ADCIRC_EqnState2(z,SP,T,zeta);
% 
% Bx = 0*zeta; By = 0*zeta; 
% for kk = 2:length(z)
%     [dpy,dpx] = gradient(pADC(:,:,kk),111e3*(lon(2)-lon(1)));
%     ii = ~isnan(dpx); 
%     %zup = 0*zeta + z(kk-1);
%     %znow = 0*zeta + z(kk); znow
%     Bx(ii) = Bx(ii) + dpx(ii)*(z(kk)-z(kk-1));
%     ii = ~isnan(dpy);
%     By(ii) = By(ii) + dpy(ii)*(z(kk)-z(kk-1));
% end

% DPX = zeros(numel(zeta),length(z));
% DPY = zeros(numel(zeta),length(z));
% for kk = 2:length(z)
%     [dpy,dpx] = gradient(pADC(:,:,kk),111e3*(lon(2)-lon(1)));
%     ii = ~isnan(dpx);
%     DPX(ii,kk) = dpx(ii);
%     ii = ~isnan(dpy);
%     DPY(ii,kk) = dpy(ii);
% end
% DPX = reshape(DPX,size(zeta,1),[],length(z));
% DPY = reshape(DPY,size(zeta,1),[],length(z));
% Bx = trapz(z,DPX,3); By = trapz(z,DPY,3); 

%jj = [1 20 28 33 40];
% Lets just take smaller sub-sample for test
%SP = SP(:,:,jj);
%T = T(:,:,jj);
%p = p(jj);