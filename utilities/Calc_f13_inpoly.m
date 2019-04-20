function obj = Calc_f13_inpoly(obj,attribute,polys,cfvals,inverse,Cf_val_on_mesh)
% obj = Calc_f13_inpoly(obj,attribute,polys,cfvals,inverse,Cf_val_on_mesh)
% Input a msh class object with bathy data and get back the specified
% attribute in a f13 structure of the obj
%
%  Inputs:   1) A msh class obj with bathy on it
%            2) The attribute indicator
%               Attributes currently supported:
%              'Cf' ('quadratic_friction_coefficient_at_sea_floor')
%              'EV' ('average_horizontal_eddy_viscosity_in_sea_water_wrt_depth')
%               etc. 
%            3) A cell-arry of polygons in which you would like to alter 
%               the attribute.
%            4) A set of attribute values that correspond 1-to-1 with the
%               cell of polygons.
%
%  Outputs: 1) msh class obj with attribute values populating the f13 struct
%
%  Author:      Keith Roberts, WP to make it for general attribute
%  Created:     April 5 2018, July 5 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    inverse = 0*cfvals;
end

if strcmpi(attribute,'Cf')
    attrname = 'quadratic_friction_coefficient_at_sea_floor';
    default_val = 0.0025; % Default Cf
elseif strcmpi(attribute,'EV')
    attrname = 'average_horizontal_eddy_viscosity_in_sea_water_wrt_depth';
    default_val = 0.05; % Default smagorinsky eddy viscosity
elseif strcmpi(attribute,'mannings_n_at_sea_floor')
    attrname = 'mannings_n_at_sea_floor';
    default_val = 0.020;
elseif strcmpi(attribute,'surface_submergence_state')
    attrname = 'surface_submergence_state';
    default_val = 0;
else
    error(['Attribute ' attribute ' not currently supported'])
end

%% Make into f13 struct
if isempty(obj.f13)
    % Add internal_tide as first entry in f13 struct
    obj.f13.AGRID = obj.title;
    obj.f13.NumOfNodes = length(obj.b);
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
if nargin < 6
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
%EOF
end



