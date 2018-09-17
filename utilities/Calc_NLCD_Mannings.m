function obj = Calc_NLCD_Mannings(obj,NLCD)
% obj = Calc_f13_inpoly(obj,NCLD)
% Input a msh class object and interpolates a NLCD land-cover database file
% onto the msh while converting to mannings values through look-up table
%
% Download:
% NLCD 2006 Land Cover (2011 Edition) (1.1Gb) from
% https://www.mrlc.gov/nlcd06_data.php
% 
%
%  Author:   WP July 19, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NCLD table
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

% The mannings name and default value
attrname = 'mannings_n_at_sea_floor';
default_val = 0.02;
dmy = msh();  dmy.p = obj.p; dmy.t = obj.t; 
% Get NCLD values on dummy msh object using nearest neighbour
obj1 = interp(dmy,NLCD,'interp','nearest','type','depth');
obj1.b(isnan(obj1.b),:)=11; % <--default value to NaN 
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



