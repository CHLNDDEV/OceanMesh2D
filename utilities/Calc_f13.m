function obj = Calc_f13(obj,attribute,varargin) 
% obj = Calc_f13(obj,attribute,varargin)
% Input a msh class object and get back the specified
% attribute in a f13 structure of the obj
%
%  Inputs:   1) A msh class obj with bathy on it
%            2) The attribute indicator
%               Attributes currently supported:
%              'Cf' ('quadratic_friction_coefficient_at_sea_floor')
%              'Ev' ('average_horizontal_eddy_viscosity_in_sea_water_wrt_depth')
%              'Mn' ('mannings_n_at_sea_floor')
%              'Ss' ('surface_submergence_state')
%              'Re' ('initial_river_elevation')
%              'Ad' ('advection_state')
%              'Sb' ('subgrid_barrier')
%              'Es' ('elemental_slope_limiter')
%
%            3) then either: 
%              'inpoly' followed by...
%            - A cell-arry of polygons in which you would like to alter 
%               the attribute.
%            - A set of attribute values that correspond 1-to-1 with the
%               cell of polygons.
%            - (optional) A set of 0 or 1's that correspond 1-to-1 with the
%               cell of polygons as to whether the in or out polygon is
%               selected
%
%            or 'assign' followed by...
%            - an array of values to assign with length the same as the 
%              number of vertices in the msh obj 
%
%  Outputs: 1) msh class obj with attribute values populating the f13 struct
%
%  Author:      Keith Roberts, WP to make it for general attribute
%  Created:     April 5 2018, July 5 2018, June 6 2019 (cleaning up)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(attribute,'Cf')
    attrname = 'quadratic_friction_coefficient_at_sea_floor';
    default_val = 0.0025; % Default Cf
elseif strcmpi(attribute,'Ev')
    attrname = 'average_horizontal_eddy_viscosity_in_sea_water_wrt_depth';
    default_val = 0.05; % Default smagorinsky eddy viscosity
elseif strcmpi(attribute,'Mn')
    attrname = 'mannings_n_at_sea_floor';
    default_val = 0.020;
elseif strcmpi(attribute,'Ss')
    attrname = 'surface_submergence_state';
    default_val = 0;
elseif strcmpi(attribute,'Re')
    attrname = 'initial_river_elevation';
    default_val = 0;
elseif strcmpi(attribute,'Ad')
    attrname = 'advection_state';
    default_val = -999;
elseif strcmpi(attribute,'Sb')
     attrname = 'subgrid_barrier';
     default_val = 99999;
elseif strcmpi(attribute,'Es')
     attrname = 'elemental_slope_limiter';
     default_val = -1e-3;
else
    error(['Attribute ' attribute ' not currently supported'])
end

if strcmpi(varargin{1},'inpoly')
    polys = varargin{2};
    cfvals = varargin{3};
    if length(varargin) < 4
        inverse = 0*cfvals;
    else
        inverse = varargin{4};
    end
elseif strcmpi(varargin{1},'assign')
     Cf_val_on_mesh = varargin{2};
else
    error(['First varargin entry (' varargin{1} ') unsupported'])
end

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
            % overwrite existing tau0
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
cf = obj.p(:,1)*0 + default_val;
if strcmpi(varargin{1},'inpoly')
    for i = 1 : length(polys)
        in = inpoly([obj.p(:,1),obj.p(:,2)],polys{i});
        if inverse(i)
            cf(~in) = cfvals(i);
        else
            cf(in) = cfvals(i);
        end
    end
else
    cf = Cf_val_on_mesh;
end
numnodes = length(find(cf ~= default_val));
obj.f13.userval.Atr(NA).usernumnodes = numnodes ;
% Print out list of nodes for each
K = find(cf ~= default_val);
obj.f13.userval.Atr(NA).Val = [K cf(K)]';

if ~isempty(obj.f15)
    % Change attribute in obj.f15
    disp(['Adding on ' attrname ' into fort.15 struct'])
    obj.f15.nwp = obj.f15.nwp + 1;
    obj.f15.AttrName(obj.f15.nwp).name = attrname;
end

%EOF
end
