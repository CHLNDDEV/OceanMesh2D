function obj = Calc_tau0(obj,varargin)
% obj = Calc_tau0(obj,varargin)
% Input a msh class object with bathy data and get back the tau0 primitive
% weighting for continuity equation in a f13 structure of the obj
%
%  Inputs:   1) A msh class obj with bathy on it
%            2) Optional name-value arguments:
%             - 'opt': -3 [default] => calculates spatially varying tau0
%                                      factors
%                     +ve => calculates the stable positive tau0 based on 
%                            the value of 'opt', which is the intended 
%                            simulation timestep
%
%             additional varargin options for 'opt' = -3:
%             - 'distance': the cutoff mean distance to neighbouring nodes
%               to switch between depth based tau0 (below) and the default
%               value, 0.03 (default distance is 2 [km])
%             - 'depth': the cutoff depth to switch between 0.005
%               and 0.02 (default is 10 [m])
%
%             additional varargin options for 'opt' = +ve:
%             - 'kappa': the value of the time weighting factor, kappa for 
%                        the future time step. Must be <= 1 and is 0.5 
%                        by default.
%             - 'sf': the safety factor applied to the stable value of
%                     tau0. Must be <= 1 and is 0.6 by default. 
%
%  Outputs: 1) msh class obj with primitive_weighting_in_continuity_equation
%              values populating the f13 struct
%
%  Author:      William Pringle                                 
%  Created:     March 14 2018                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

attrname = 'primitive_weighting_in_continuity_equation';

%% Test optional arguments
% default
MinDepth = 10; % m 
Distance = 2; % km
opt = -3; % calculate tau0 using -3 option
mm  = 2/3; % consistent mass matrix
kappa = 0.4; % GWCE weight
sf = 0.6; % suggested safety factor
if ~isempty(varargin)
    names = {'depth','distance','opt','kappa','sf'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                MinDepth = varargin{ind*2};
            elseif ii == 2
                Distance = varargin{ind*2};
            elseif ii == 3
                opt = varargin{ind*2};
            elseif ii == 4
                kappa = varargin{ind*2};
            elseif ii == 5
                sf = varargin{ind*2};
            end
        end    
    end
end

if opt > 0
   dt = opt;
   obj.f15.tau0 = sf*4*(2-mm)*(3*kappa-1)/dt;
   obj.f15.a00b00c00 = [kappa kappa 1-2*kappa];
   obj.f15.im = 511111;
   disp(['Computed constant value of tau0 = ' num2str(obj.f15.tau0) ' based on...'])
   disp(['mm = ' num2str(mm)])
   disp(['kappa = ' num2str(kappa)])
   disp(['sf = ' num2str(sf)])
   disp(['dt = ' num2str(dt)])
   return;
end

if isempty(obj.b) || all(obj.b == 0)
    error('No bathymetry data in grid to calculate the tau0 coefficients')
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



