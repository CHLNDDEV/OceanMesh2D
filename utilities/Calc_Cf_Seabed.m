function obj = Calc_Cf_Seabed(obj,datatables,uf,depth_cutoff)
% obj = Calc_Cf_Seabed(obj,datatables,uf,depth_cutoff)
%
% Put Cf in f13 struct calculated using usSEABED TXT datatables
% You can download these files from:
% https://coastalmap.marine.usgs.gov/national/usseabed/data.html
% I recommend selecting the clc (calculated data - download the txt file)
% 
% Formula used is:
% 
% Cf = [kappa/ln(0.5*h/z0)]^2
% where h is water depth, kappa = 0.4, and z0 is calculated using the van
% Rijn (2007) formula. Procedure to get z0 is detailed in my publication:
%
% William J. Pringle, et al. Finite-Element barotropic model for the Indian
% and Western Pacific Oceans: Tidal model-data comparisons and sensitivities,
% Ocean Modelling, Volume 129, 2018 (13-38), doi: 10.1016/j.ocemod.2018.07.003.
%
% Inputs:
%       obj: msh class obj
%       datatables: list of txt filenames of the usSEABED data
%       uf: vector of the flow velocities to use - must be same length as the
%           vertices in the msh class obj
%       depth_cutoff: a 2-tuple for lower and upper limit on the depth, 
%           respectively. The upper limit is where we switch to default Cf 
%           of 0.0025
%
% Output: msh class obj with the f13 struct populated for the 
%         quadratic_friction_coefficient_at_sea_floor attribute
%
% Creation: 8/26-27/2019 William Pringle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The attribute name and default value we are using
attrname = 'quadratic_friction_coefficient_at_sea_floor';
default_val = 0.0025; % Default Cf

%% Check inputs
if isempty(obj.b) || all(obj.b == 0)
   error('Depths are required to be on the msh') 
end

% set depth into h vector
h = obj.b;
h = max(h,depth_cutoff(1)); % make sure at least depth_cutoff(1) deep

if length(uf) ~= 1 && length(uf) ~= length(h)
   error('Input velocity vector uf must have same length as mesh vertices')
end 

%% Read and Process the usSeabed data
Loc = []; Phi = []; CrShSt = []; RH = [];
for tn = datatables
    T = readtable(tn);
    Loc = [Loc; T.Longitude T.Latitude];
    Phi =  [Phi; T.Grainsize];
    if sum(strcmp('Lgcrshst',fieldnames(T))) > 0
        CrShSt = [CrShSt; T.Lgcrshst];
    elseif sum(strcmp('LgCrShSt',fieldnames(T))) > 0
        CrShSt = [CrShSt; T.LgCrShSt];
    else
        error('No crit chear stress field in datatable')
    end
    RH = [RH; T.Roughness];
end

% Get rid of null data and extract the vertical roughness height only
VRLoc = Loc(~strcmp(RH,'-99') & ~cellfun('isempty',RH),:);
RH = RH(~strcmp(RH,'-99') & ~cellfun('isempty',RH));
VR = 0*VRLoc(:,1);
for ii = 1:length(RH)
    RHchar = RH{ii};
    VR(ii) = str2double(RHchar(1:strfind(RHchar,':')-1));
end

% Make the scatteredInterpolant for Phi & Shear stress
% note that null data is value equal to -99
Fp = scatteredInterpolant(Loc(Phi > -99,:),Phi(Phi > -99),'nearest','none');
Fc = scatteredInterpolant(Loc(CrShSt > -99,:),CrShSt(CrShSt > -99),'nearest','none');
Fv = scatteredInterpolant(VRLoc,VR,'nearest','none');

%% Set the sediment and flow properties
dgrav = 2e-3;
dsand = 6.2e-5;
dsilt = 3.2e-5;
D0 = 1e-3; % the reference diameter
d50 = D0*2.^(-Fp(obj.p)); % get the median sediment diameter in meters
Ucrit2 = 10.^Fc(obj.p); % get the critical shear stress in kPa 
VR = 2.^Fv(obj.p)/100; % vertical roughness height in meters

% Calculate mobility parameter from the knowledge of critical shear stress
% (avoids having to guess the density)
mob = uf.^2./Ucrit2;

%% Now Follow Appendix A from Pringle et al., 2018 to get ks roughness
% get the ripple and mega-ripple roughness factors
fsc = min(1,(0.25*dgrav./d50).^(1.5));
ffs = min(1,d50./(1.5*dsand));

% get the individual ripple, mega-rippled and dune roughnesses
ksr = fsc.*d50.*(85-65*tanh(0.0015*(mob-150)));
ksm = max(min(0.02,200*d50),2e-5*ffs.*h.*(1-exp(-0.05*mob)).*(500-mob));
ksd = max(0,8e-5*ffs.*h.*(1-exp(-0.02*mob)).*(600-mob));

% get the total equivalent sand roughness
ks = max(20*dsilt,sqrt(ksr.^2 + ksm.^2 + ksd.^2));

%% Convert to Cf
% Max of the ks or the vertical roughness height
ks = max(ks,VR);
z0 = ks/30;
kappa = 0.4;
Cf = (kappa./log(0.5*h./z0)).^2;
% just set deep areas to default roughness
Cf(h > depth_cutoff(2) | isnan(d50)) = default_val;

%% Make into f13 struct
if isempty(obj.f13)
    % Add this as first entry in f13 struct
    obj.f13.AGRID = obj.title;
    obj.f13.NumOfNodes = length(obj.p);
    obj.f13.nAttr = 1;
    NA = 1;
else
    broken = 0;
    for NA = 1:obj.f13.nAttr
        if strcmp(attrname,obj.f13.defval.Atr(NA).AttrName)
            broken = 1;
            % overwrite existing Cf
            break
        end
    end
    if ~broken
        % add internal_tide to list
        obj.f13.nAttr = obj.f13.nAttr + 1;
        NA = obj.f13.nAttr;
    end
end

% Default Values
obj.f13.defval.Atr(NA).AttrName = attrname;
% We can just put in the options here
obj.f13.defval.Atr(NA).Unit = 'unitless';
valpernode = 1;
obj.f13.defval.Atr(NA).ValuesPerNode = valpernode ;
obj.f13.defval.Atr(NA).Val = default_val ;

% User Values
obj.f13.userval.Atr(NA).AttrName = attrname ;
numnodes = length(find(Cf ~= default_val));
obj.f13.userval.Atr(NA).usernumnodes = numnodes ;
% Print out list of nodes for each
K = find(Cf ~= default_val);
obj.f13.userval.Atr(NA).Val = [K Cf(K)]';

end
