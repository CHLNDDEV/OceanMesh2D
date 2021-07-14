function obj = Calc_Mannings_Landcover(obj,data,type,varargin)
% obj = Calc_Mannings_Landcover(obj,data,type,varargin)
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
% varargin: Accepts the same options as for msh.interp to control how 
%           data is interpolated; see 'help msh.interp'
% 
%  Author:            KR & WP July 19, 2018
%  Update for CCAP:   WP May 5, 2021
%  Update for storing incremental interpolation results: KJ & OK July 9, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(type)
   type = 'nlcd'; 
end

if strcmp(type,'nlcd')
    disp('Info: Using NLCD table')
    % NLCD table
    load nlcd
    nlcd_mannings(end+1)=nlcd_mannings(11); % <--make NaN = default value
    varargin{end+1}='lut';
    varargin{end+1}=nlcd_mannings;
elseif strcmp(type,'ccap')
    disp('Info: Using CCAP table')
    % CCAP table
    load ccap 
    ccap_mannings(end+1)=ccap_mannings(21); % <-- make NaN = default value
    varargin{end+1}='lut';
    varargin{end+1}=ccap_mannings;
else
    error('Land-cover database not supported')
end
% The mannings name and default value
attrname = 'mannings_n_at_sea_floor';
default_val = varargin{end}(end); %end of the lut is default
dmy = msh();  dmy.p = obj.p; dmy.t = obj.t; 
% Convert to Mannings and interpolate how the user wants
dmy = interp(dmy,data,varargin{:});
Man = dmy.b';
Man( isnan(Man) | Man==0 ) = default_val; %points outside of lut can still be NaN
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
if isfield(obj.f13.userval.Atr(NA),'Val')
    obj.f13.userval.Atr(NA).Val = [obj.f13.userval.Atr(NA).Val'; K' , Man(K)']';
else
    obj.f13.userval.Atr(NA).Val = [K ; Man(K)];
end
%EOF
end


