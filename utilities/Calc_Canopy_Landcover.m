function obj = Calc_Canopy_Landcover(obj,data,type,varargin)
% obj = Calc_Canopy_Landcover(obj,data,type,varargin)
% Input a msh class object and interpolates a land-cover database file
% (netcdf format) onto the msh while converting to surface canopy
% coefficient through a look-up table
% 
% type:
%  - 'nlcd': 
%    NLCD 2006 Land Cover (2011 Edition) (1.1Gb) from
%    https://www.mrlc.gov/nlcd06_data.php
%  - 'ccap':
%    CCAP NOAA land cover:
%    https://coast.noaa.gov/data/digitalcoast/pdf/ccap-class-scheme-regional.pdf
%
% varargin: Accepts the same options as for msh.interp to control how 
%           data is interpolated; see 'help msh.interp'
%
% Usage:
%   m = Calc_Canopy_Landcover(m,landuse_data,'ccap','interp','nearest');
%   m = Calc_Canopy_Landcover(m,landuse_data,'nlcd','interp','CA','N',2);
% 
%  Author:            WP Jun, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(type)
   type = 'nlcd'; 
end

if strcmp(type,'nlcd')
    disp('Info: Using NLCD table')
    % NLCD table
    load nlcd
    nlcd_canopy(end+1)=nlcd_canopy(11); % <--make NaN = default value
    varargin{end+1}='lut';
    varargin{end+1}=nlcd_canopy;
elseif strcmp(type,'ccap')
    disp('Info: Using CCAP table')
    % CCAP table
    load ccap 
    ccap_canopy(end+1)=ccap_canopy(21); % <-- make NaN = default value
    varargin{end+1}='lut';
    varargin{end+1}=ccap_canopy;
else
    error('Land-cover database not supported')
end
% The canopy name and default value
attrname = 'surface_canopy_coefficient';
default_val = 1; % (no canopy)
dmy = msh();  dmy.p = obj.p; dmy.t = obj.t; 
% Convert to surface canopy and interpolate how the user wants
dmy = interp(dmy,data,varargin{:});
Can = dmy.b';
Can( isnan(Can) | Can==0 ) = default_val; %points outside of lut can still be NaN
%% Make into f13 struct
if isempty(obj.f13)
    % Add add canopy as first entry in f13 struct
    obj.f13.AGRID = obj.title;
    obj.f13.NumOfNodes = length(obj.p);
    obj.f13.nAttr = 1;
    NA = 1;
else
    broken = 0;
    for NA = 1:obj.f13.nAttr
        if strcmp(attrname,obj.f13.defval.Atr(NA).AttrName)
            broken = 1;
            % overwrite existing canopy
            break
        end
    end
    if ~broken
        % add canopy to list
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
numnodes = length(find(Can ~= default_val));
obj.f13.userval.Atr(NA).usernumnodes = numnodes ;
% Print out list of nodes for each
K = find(Can ~= default_val);
obj.f13.userval.Atr(NA).Val = [K ; Can(K)];
%EOF
end
