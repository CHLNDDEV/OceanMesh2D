function obj = tide_fac(obj,t_s,t_e,incnstit)
% obj = tide_fac(obj,t_s,t_e,incnstit)
% Input a msh class obj, start and time simulation times, and constituent
% names and return populated f15 tidal potential options in the msh class
% obj                
%                                                                       
% Created by William Pringle March 15 2018, 
% Copied functions from U_tide to calculate the tide factors, available at: 
% https://www.mathworks.com/matlabcentral/fileexchange/
% 46523--utide--unified-tidal-analysis-and-prediction-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the times and lat entries
ts = datenum(t_s);
te = datenum(t_e);

tref = 0.5*(ts + te);

% make a fake t with 1000 intervals
t = linspace(ts,te,1000)';

lat = mean(obj.p(:,2));

%% Allow for an entry just to be major8 to get the major eight constituents
if strcmp(incnstit{1},'major8')
   incnstit = {'M2','S2','N2','K2','K1','O1','Q1','P1'};
   disp('Using the major eight harmonic tidal constituents')
   disp(incnstit)
end

%% Now check and get the indices of the incnstits (this part came from U_tide)
load('tide_fac_constants.mat','const');
cnstit.NR.lind = nan*ones(size(incnstit,1),1);
for j = 1:length(incnstit)
    lind1 = strmatch(incnstit{j},const.name);
    if isempty(lind1)
        error(['tide_fac: unrecognized non-reference constituent: '...
               incnstit{j} '.']);
    else
        cnstit.NR.lind(j) = lind1;
    end
end
[~,seq] = sort(const.freq(cnstit.NR.lind));
cnstit.NR.lind = cnstit.NR.lind(seq);
if isempty(cnstit.NR.lind)
    error('tide_fac: no constituents specified');
end
cnstit.NR.frq = const.freq(cnstit.NR.lind);
cnstit.NR.name = cellstr(const.name(cnstit.NR.lind,:));
% Doodson amp multiplied by Doodson constant/g:  
% (See : https://www.whoi.edu/fileserver.do?id=21351&pt=10&p=17272)
cnstit.NR.amp  = abs(const.doodsonamp(cnstit.NR.lind))*0.2675; 
cnstit.NR.red  = const.earthreduc(cnstit.NR.lind); 

%% Now get the F, U, and V (using the U_tide function ut_FUV)
[F,U,V] = ut_FUV(t,tref,cnstit.NR.lind,lat,zeros(1,4));

%% Final factors
F = mean(F)'; % This is the average nodal factor
phs = (mean(U)+V(1,:))'*360;  % This is the average nodal correction + 
                              % astronomical argument at beginning of
                              % simulation
                              
% Make sure phase between 0 and 360 for aesthetic purposes
while any(phs < 0)
    phs(phs<0) = phs(phs<0) + 360;
end

obj.f15.ntif = length(cnstit.NR.name);
for k = 1:obj.f15.ntif
    obj.f15.tipotag(k).name   = cnstit.NR.name{k}; % name in frq order
    obj.f15.tipotag(k).val(1) = cnstit.NR.amp(k); % potential amplitude of species
    obj.f15.tipotag(k).val(2) = cnstit.NR.frq(k)*2*pi/3600; % frq in rad/s format
    obj.f15.tipotag(k).val(3) = cnstit.NR.red(k); % earth rigidity reduction factor
    obj.f15.tipotag(k).val(4) = F(k);   % average nodal factor
    obj.f15.tipotag(k).val(5) = phs(k); % average nodal correction + 
                         % astronomical argument at beginning of simulation
end
%EOF
end

% Stuff below has been copied and slightly adjusted from U_tide 
% (which was copied from T_tide)
function [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs)
% UT_FUV()
% compute nodal/satellite correction factors and astronomical argument
% inputs
%   t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%   lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
% output
%   F = real nodsat correction to amplitude [unitless] (nt x nc)
%   U = nodsat correction to phase [cycles] (nt x nc)
%   V = astronomical argument [cycles] (nt x nc)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)

nt = length(t);
nc = length(lind);
%% nodsat
if ngflgs(2) % none
    F = ones(nt,nc);
    U = zeros(nt,nc);
else
    if ngflgs(1) % linearized times
        tt = tref;
    else         % exact times
        tt = t;
    end
    ntt = length(tt);
    load('tide_fac_constants.mat');
    [astro,~]=ut_astron(tt');
    if abs(lat) < 0.5 % to make sure no close to zero values
        lat = sign(lat)*0.5; 
        if sign(lat) == 0
            lat = 0.5;
        end
    end
    slat=sin(pi*lat/180);
    rr=sat.amprat;
    j=find(sat.ilatfac==1);
    rr(j)=rr(j).*0.36309.*(1.0-5.0.*slat.*slat)./slat;
    j=find(sat.ilatfac==2);
    rr(j)=rr(j).*2.59808.*slat; 
    uu=rem( sat.deldood*astro(4:6,:)+sat.phcorr(:,ones(1,ntt)), 1);
    nfreq=length(const.isat); %#ok
    mat = rr(:,ones(1,ntt)).*exp(1i*2*pi*uu);
    F = ones(nfreq,ntt);
    ind = unique(sat.iconst);
    for i = 1:length(ind)
        F(ind(i),:) = 1+sum(mat(sat.iconst==ind(i),:),1);
    end
    U = imag(log(F))/(2*pi); % faster than angle(F)
    F=abs(F);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        exp2 = abs(exp1);
        F(k,:)=prod(F(j,:).^exp2(:,ones(ntt,1)),1);
        U(k,:)=sum(U(j,:).*exp1(:,ones(ntt,1)),1);
    end
    F=F(lind,:)';
    U=U(lind,:)';
    if ngflgs(1) % nodal/satellite with linearized times
        F = F(ones(nt,1),:);
        U = U(ones(nt,1),:);
    end
end
%% gwch (astron arg)
if ngflgs(4) % none (raw phase lags not greenwich phase lags)
    if ~exist('const','var')
        load('tide_fac_constants.mat','const');
    end
    [~,ader] = ut_astron(tref);
    ii=isfinite(const.ishallow); 
    const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
    for k=find(ii)'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
    end
    V = 24*(t-tref)*const.freq(lind)';
else 
    if ngflgs(3)  % linearized times
        tt = tref;
    else 
        tt = t;   % exact times
    end
    ntt = length(tt);
    if exist('astro','var')
        if ~isequal(size(astro,2),ntt)
            [astro,~]=ut_astron(tt');
        end        
    else
        [astro,~]=ut_astron(tt');
    end
    if ~exist('const','var')
        load('tide_fac_constants.mat');
    end
    V=rem( const.doodson*astro+const.semi(:,ones(1,ntt)), 1);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        V(k,:) = sum(V(j,:).*exp1(:,ones(ntt,1)),1);
    end
    V=V(lind,:)';
    if ngflgs(3)    % linearized times
        [~,ader] = ut_astron(tref);
        ii=isfinite(const.ishallow);
        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
        for k=find(ii)'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
                shallow.coef(ik) );
        end
        V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
    end
end
%EOF
end

function [astro,ader] = ut_astron(jd)
% UT_ASTRON()
% calculate astronomical constants
% input
%   jd = time [datenum UTC] (1 x nt)
% outputs
%   astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
%   ader = matrix of derivatives of astro [cycles/day] (6 x nt)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (copy of t_astron.m from t_tide, Pawlowicz et al 2002)

d=jd(:)'-datenum(1899,12,31,12,0,0);
D=d/10000;
args=[ones(size(jd));
      d;
      D.*D;
      D.^3];
sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];
astro=rem( [sc;hc;pc;npc;ppc]*args./360.0 ,1);
tau=rem(jd(:)',1)+astro(2,:)-astro(1,:);
astro=[tau;astro];
dargs=[zeros(size(jd));
       ones(size(jd));
       2.0e-4.*D;
       3.0e-4.*D.*D];
ader=[sc;hc;pc;npc;ppc]*dargs./360.0;
dtau=1.0+ader(2,:)-ader(1,:);
ader=[dtau;ader];
end