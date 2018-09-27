function obj = Calc_tau0(obj,varargin)
% obj = Calc_tau0(obj,varargin)
% Input a msh class object with bathy data and get back the tau0 primitive
% weighting for continuity equation in a f13 structure of the obj
%
%  Inputs:   1) A msh class obj with bathy on it
%            2) Optional name-value arguments:
%               'distance': the cutoff mean distance to neighbouring nodes
%               to switch between depth based tau0 (below) and the default
%               value, 0.03 (default distance is 2 [km])
%               'depth': the cutoff depth to switch between 0.005
%               and 0.02 (default is 10 [m])
%
%  Outputs: 1) msh class obj with primitive_weighting_in_continuity_equation
%              values populating the f13 struct
%
%  Author:      William Pringle                                 
%  Created:     March 14 2018                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

attrname = 'primitive_weighting_in_continuity_equation';

if isempty(obj.b)
    error('No bathymetry data in grid to calculate the tau0 coefficients')
end

%% Test optional arguments
% default
MinDepth = 10; % m 
Distance = 2; % km
if ~isempty(varargin)
    names = {'depth','distance'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                MinDepth = varargin{ind*2};
            elseif ii == 2
                Distance = varargin{ind*2};
            end
        end    
    end
end

%% Get all the unique bar edges and distances
all_bars = [obj.t(:,[1,2]); obj.t(:,[1,3]); obj.t(:,[2,3])];   
bars = unique(sort(all_bars,2),'rows'); 
nodepairs = reshape(bars',[],1);
distances = m_lldist(obj.p(nodepairs,1),obj.p(nodepairs,2));
distances = distances(1:2:end);

%% Sort the bars, distances, and work out the numbers of connections 
%% for each vertex.
bars = [bars; fliplr(bars)];
distances = repmat(distances,2,1);
[bars,I] = sortrows(bars,1); distances = distances(I);
edges = unique(bars);
N = histc(bars(:,1),edges);

%% Loop through and take the average distance between connecting vertices
ns = 1; Md = zeros(size(obj.b)); 
for ii = 1:length(obj.b)
    ne = ns + N(ii) - 1; 
    Md(ii) = mean(distances(ns:ne));
    ns = ne + 1;
end

%% Calculate tau0 based on conditions
default_val = 0.03;
tau0 = default_val*ones(size(obj.b));
tau0(Md > Distance & obj.b <= MinDepth) = 0.02;
tau0(Md > Distance & obj.b > MinDepth) = 0.005;

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
numnodes = length(find(tau0 < default_val));
obj.f13.userval.Atr(NA).usernumnodes = numnodes ;
% Print out list of nodes for each
K = find(tau0 < default_val);
obj.f13.userval.Atr(NA).Val = [K tau0(K)]';

if ~isempty(obj.f15)
    % Change attribute in obj.f15
    disp('Changing tau0 parameter to -3')
    disp('Adding tau0 attribute name in fort.15 struct')
    obj.f15.tau0 = -3;
    obj.f15.nwp = obj.f15.nwp + 1;
    obj.f15.AttrName(obj.f15.nwp).name = attrname;
end
%EOF
end 



