function obj = Calc_Sponge(obj,W,generator,varargin)
% obj = Calc_Sponge(obj,W,generator,varargin)
% Adds the sponge information into the f13 struct in the msh class obj.
% Must specify the width (W) of the sponge and the sponge layer type:
% = 1 -> generating-absorbing with tidal solutions (need fort.53001/54001)
% = 0 -> fully absorbing sponge 
% = -1 -> generating-absorbing with arbitrary solutions (need fort.2001)
% Other inputs are optional.
%
% Default for other inputs are:
% F = 20, dis_ratio = 0.5 and spngtype = 'poly', 
% i.e. that a polynomial (2nd order) type function on sigma will be used so 
% that half of the way into the sponge (dis_ratio) the signal will be 
% attenuated 20 fold (F)

if isempty(obj.b)
   error('This function requires bathymetry on the msh object (use interp)') 
end

%% Test optional arguments
% default
g         = 9.81;  % gravity
F         = 20;    % 20 fold decrease at dis_ratio
dis_ratio = 0.5;   % ratio of sponge width where F decrease occurs
spngtype  = 'poly'; alpha = 2; % second order polynomial type function 
if ~isempty(varargin)
    names = {'F','dis_ratio','spngtype'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                F = varargin{ind*2}; 
            elseif ii == 2
                dis_ratio = varargin{ind*2};
            elseif ii == 3
                spngtype = varargin{ind*2};
            end
        end    
    end
end   
    
%% Get the nodes for each sponge zone
sponge = get_latlon_for_sponge_zone(obj,W);

%% Make the sigma based on spongetype and coefficients
sigma = []; idspg_node = [];
for op = 1:obj.op.nope
    % need to get metre equivalents of L & d (roughly)
    W = sponge(op).W*111e3; d = sponge(op).d*111e3;
    % spongetype polynomial or hyperbola
    if strcmp(spngtype,'poly')
        sigma_m = -sqrt(g*sponge(op).H)*(alpha+1)*log(1/F)/...
                        (W*dis_ratio^(alpha+1));  
        sigma_n = sigma_m*(d/W).^alpha;
    elseif  strcmp(spngtype,'hyper')
        sigma_m = -log(1/F)/(-log(1-dis_ratio)-dis_ratio);
        sigma_n = sigma_m*min(1,sqrt(g*sponge(op).H)*d./(W*(W-d)));
    end
    idspg_node = [idspg_node; sponge(op).idx];
    sigma = [sigma; sigma_n];
end

% The following is needed when we have multiple open boundaries in order to
% take the maximum value 
unique_nodes = unique(idspg_node);
if length(unique_nodes) < length(idspg_node)
    [idspg_node,IA] = sort(idspg_node);
    sigma = sigma(IA);
    sigma_unique = 0*unique_nodes;
    nn = 1;
    for ii = 1:length(unique_nodes)
        for jj = 1:obj.op.nope
            if idspg_node(nn) == unique_nodes(ii)
                sigma_unique(ii) = max(sigma(nn),sigma_unique(ii));
                nn = nn + 1;
                if nn > length(idspg_node)
                    break;
                end
            end
        end
    end
    idspg_node = unique_nodes;
    sigma = sigma_unique;
end
  
%% Set f13 structure in the msh object class
if isempty(obj.f13)
    % Add sponge_generator_layer as first entry in f13 struct
    obj.f13.AGRID = obj.title;
    obj.f13.NumOfNodes = length(obj.b); 
    obj.f13.nAttr = 1;
    natb = 1;
else
    broken = 0;
    for natb = 1:obj.f13.nAttr
        if strcmp('sponge_generator_layer',obj.f13.defval.Atr(natb).AttrName)
            broken = 1;
            % overwrite existing internal_tide
            break
        end
    end
    if ~broken 
        % add sponge_generator_layer to list
        obj.f13.nAttr = obj.f13.nAttr + 1;
        natb = obj.f13.nAttr;
    end
end
attrname = 'sponge_generator_layer' ; 

% User-defined input
obj.f13.userval.Atr(natb).AttrName = attrname ; 
obj.f13.userval.Atr(natb).usernumnodes = length(idspg_node) ;
obj.f13.userval.Atr(natb).Val = [ idspg_node'; sigma'; ...
                                  (sigma*0 + generator)' ] ;                          
%                   
% Default input
obj.f13.defval.Atr(natb).AttrName =  attrname ; 
obj.f13.defval.Atr(natb).Unit = 'unitless' ;
obj.f13.defval.Atr(natb).ValuesPerNode = 2 ;
obj.f13.defval.Atr(natb).Val = [0.0 0] ;

if ~isempty(obj.f15)
    % Change attribute in obj.f15
    disp('Adding on sponge attribute name in fort.15 struct')
    obj.f15.nwp = obj.f15.nwp + 1;
    obj.f15.AttrName(obj.f15.nwp).name = attrname;
end

%EOF
end

function sponge = get_latlon_for_sponge_zone(obj,W)
% sponge = get_latlon_for_sponge_zone(obj,W)
% This function reads the open boundaries of the mesh and gets the nodes in
% the sponge layer as specified by the width, W
% Inputs:
% W      = width in degrees
%          or  
%        = [period, frac]; 
% where period of a wave [s] (e.g. M2), and frac is fraction of the
% wavelength to set as sponge width W.

g = 9.81; % gravity

% Now loop over boundaries and get the sponge layer points
sponge = struct([]);
for op = 1:obj.op.nope
    nodes = obj.op.nbdv(1:obj.op.nvdll(op),op);
    % get mean depth along boundary
    H = mean(obj.b(nodes));
    if length(W) > 1
        % Get the length of the sponge (in metres)
        W1 = W(2)*W(1)*sqrt(g*H); 
        % Change length to degrees (roughly)
        W1 = W1/111e3;
    else
        W1 = W;
    end
    % Find all points W distance from the open boundary
    % Get the nodes in the trivial box from OB
    xmax = max(obj.p(nodes,1)) + W1; xmin = min(obj.p(nodes,1)) - W1;
    ymax = max(obj.p(nodes,2)) + W1; ymin = min(obj.p(nodes,2)) - W1;
    I = find(obj.p(:,2) <= ymax & obj.p(:,2) >= ymin & ...
             obj.p(:,1) <= xmax & obj.p(:,1) >= xmin);   
    % Do dsegment on all the I nodes to get distance to the OB
    [~,d] = ourKNNsearch(obj.p(nodes,:)',obj.p(I,:)',1);
    idx = I(d < W1);
    sponge(op).idx = idx;
    sponge(op).d   = W1 - d( d < W1);
    sponge(op).H   = H;
    sponge(op).W   = W1;
end

%EOF
end
