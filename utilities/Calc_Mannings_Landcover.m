function obj = Calc_Mannings_Landcover(obj,data,type)
% obj = Calc_Mannings_Landcover(obj,data,type)
% Input a msh class object and interpolates a land-cover database file
% (netcdf format) onto the msh while converting to mannings values through
% a look-up table
% 
% type:
%  - 'nlcd': 
%    NLCD 2006 Land Cover (2011 Edition) (1.1Gb) from
%    https://www.mrlc.gov/nlcd06_data.php
%  - 'ccap':
%    CCAP NOAA land cover:
%    https://coast.noaa.gov/data/digitalcoast/pdf/ccap-class-scheme-regional.pdf
%
%  Author:            WP July 19, 2018
%  Update for CCAP:   WP May 5, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(type)
   type = 'nlcd'; 
end

if strcmp(type,'nlcd')
    disp('Notice: Using NLCD table')
    % NLCD table
    def_index = 11;
    nlcd_class(11) = 0.02;    %Open Water
    nlcd_class(12) = 0.010;   %Perennial Ice/Snow
    nlcd_class(21) = 0.020;   %Developed - Open Space
    nlcd_class(22) = 0.050;   %Developed - Low Intensity
    nlcd_class(23) = 0.100;   %Developed - Medium Intensity
    nlcd_class(24) = 0.150;   %Developed - High Intensity
    nlcd_class(31) = 0.090;   %Barren Land (Rock/Sand/Clay)
    nlcd_class(32) = 0.040;   %Unconsolidated Shore
    nlcd_class(41) = 0.100;   %Deciduous Forest
    nlcd_class(42) = 0.110;   %Evergreen Forest
    nlcd_class(43) = 0.100;   %Mixed Forest
    nlcd_class(51) = 0.040;   %Dwarf Scrub
    nlcd_class(52) = 0.050;   %Shrub/Scrub
    nlcd_class(71) = 0.034;   %Grassland/Herbaceous
    nlcd_class(72) = 0.030;   %Sedge/Herbaceous
    nlcd_class(73) = 0.027;   %Lichens
    nlcd_class(74) = 0.025;   %Moss
    nlcd_class(81) = 0.033;   %Pasture/Hay
    nlcd_class(82) = 0.037;   %Cultivated Crops
    nlcd_class(90) = 0.100;   %Woody Wetlands
    nlcd_class(91) = 0.100;   %Palustrine Forested Wetland
    nlcd_class(92) = 0.048;   %Palustrine Scrub/Shrib Wetland
    nlcd_class(93) = 0.100;   %Estuarine Forested Wetland
    nlcd_class(94) = 0.048;   %Estuarine Scrub/Shrub Wetland
    nlcd_class(95) = 0.045;   %Emergent Herbaceous Wetlands
    nlcd_class(96) = 0.045;   %Palustrine Emergent Wetland (Persistant)
    nlcd_class(97) = 0.045;   %Estuarine Emergent Wetland
    nlcd_class(98) = 0.015;   %Palustrine Aquatic Bed
elseif strcmp(type,'ccap')
    disp('Notice: Using CCAP table')
    % CCAP table
    def_index = 21;
    nlcd_class([1 21]) = 0.02;   %Open Water
    nlcd_class(5) = 0.020;    %Developed - Open Space
    nlcd_class(4) = 0.050;    %Developed - Low Intensity
    nlcd_class(3) = 0.100;    %Developed - Medium Intensity
    nlcd_class(2) = 0.150;    %Developed - High Intensity
    nlcd_class(6) = 0.037;    %Cultivated Crops
    nlcd_class(7) = 0.033;    %Pasture/Hay
    nlcd_class(12) = 0.050;   %Shrub/Scrub
    nlcd_class(8) = 0.034;    %Grassland/Herbaceous
    nlcd_class(9) = 0.100;    %Deciduous Forest
    nlcd_class(10) = 0.110;   %Evergreen Forest
    nlcd_class(11) = 0.100;   %Mixed Forest
    nlcd_class(20) = 0.090;   %Barren Land (Rock/Sand/Clay)
    nlcd_class(24) = 0.090;   %Tundra
    nlcd_class(25) = 0.010;   %Perennial Ice/Snow
    nlcd_class(13) = 0.100;   %Palustrine Forested Wetland
    nlcd_class(14) = 0.048;   %Palustrine Scrub/Shrub Wetland
    nlcd_class(15) = 0.045;   %Palustrine Emergent Wetland (Persistant)
    nlcd_class(22) = 0.015;   %Palustrine Aquatic Bed
    nlcd_class(16) = 0.100;   %Estuarine Forested Wetland
    nlcd_class(17) = 0.048;   %Estuarine Scrub/Shrub Wetland
    nlcd_class(18) = 0.045;   %Estuarine Emergent Wetland
    nlcd_class(19) = 0.040;   %Unconsolidated Shore
    nlcd_class(23) = 0.015;   %Estuarine Aquatic Bed
end
% The mannings name and default value
attrname = 'mannings_n_at_sea_floor';
default_val = 0.02;
dmy = msh();  dmy.p = obj.p; dmy.t = obj.t; 
% Get NCLD values on dummy msh object using nearest neighbour
obj1 = interp(dmy,data,'interp','nearest','type','depth');
obj1.b(isnan(obj1.b) | obj1.b == 0) = def_index; % <--default value to NaN 
% Convert to Mannings
Man = nlcd_class(abs(obj1.b));

%% Make into f13 struct
if isempty(obj.f13)
    % Add add mannings as first entry in f13 struct
    obj.f13.AGRID = obj.title;
    obj.f13.NumOfNodes = length(obj.p);
    obj.f13.nAttr = 1;
    NA = 1;
else
    broken = 0;
    for NA = 1:obj.f13.nAttr
        if strcmp(attrname,obj.f13.defval.Atr(NA).AttrName)
            broken = 1;
            % overwrite existing tau0
            break
        end
    end
    if ~broken
        % add mannings to list
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
numnodes = length(find(Man ~= default_val));
obj.f13.userval.Atr(NA).usernumnodes = numnodes ;
% Print out list of nodes for each
K = find(Man ~= default_val);
obj.f13.userval.Atr(NA).Val = [K ; Man(K)];
%EOF
end



