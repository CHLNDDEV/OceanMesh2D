function obj = Calc_IT_Fric(obj,N_filename,varargin)
% obj = Calc_IT_Fric(obj,N_filename,varargin)
% Input a msh class object with bathy and slope data, get the values of N
% over the depth based on N_filename (if is not empty) and calculate the
% internal tide friction based on the input parameters   
% 
%  Inputs:      1) A msh class obj with slopes and bathy on it
%               2) Required input (unless you are William Pringle): 
%               N_filename: filename and directory of the buoyancy
%               frequency .mat file
%               Warning: An empty filename will work in which it will just 
%               output the slopes*Cit and ADCIRC will expect instead some 
%               buoyancy frequency information given to it inside the code
%               3) Optional name-value arguments:
%               'type': 'local_dir' (default), 'local_sca', or 'nonlocal'
%               'crit': 'Nb' (default), 'Nmw'`
%               'Cit' : a multiplication factor (2.75 is default value 
%               based on Buijsman et al. (2016) and Pringle et al. (2018ab)
%               'cutoff_depth': the depth which to cut off the
%               internal_tide_friction below (100 [m] is the default)
%
%  Outputs:     1) msh class obj with internal_tide_friction values
%                  populating the f13 struct
%
%  Author:      William Pringle                                 
%  Created:     March 13 2018, Updated June 15 2017                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(obj.b) || isempty(obj.bx) || isempty(obj.by)
   error(['Must put bathymetry and bathymetric slopes on the msh ' ...
          'object (use GridData)']) 
end
% Some constants
% Choose projection type (can be any, not restricted to evenly-gridded data)
proj = 'Mercator';
% Coriolis coefficient
psi = 2*7.29212d-5;
% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency

%% Test optional arguments
% default
type = 'local_dir';
crit = 'Nb';
C_it = 2.75;
MinDepth = 100; % m 
if ~isempty(varargin)
    names = {'type','cutoff_depth','Cit','crit'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                type = varargin{ind*2}; 
            elseif ii == 2
                MinDepth = varargin{ind*2};
            elseif ii == 3
                C_it = varargin{ind*2};
            elseif ii == 4 
                crit = varargin{ind*2};
            end
        end    
    end
end   

% Coriolis
f = psi*sind(obj.p(:,2));           
    
%% Load the constant contours of N values and compute Nb and Nmean
if ~isempty(N_filename)
    load(N_filename);  
    [Nb,Nm,Nmw] = Gridded_to_Mesh_SeaBed_DepthAveraged(obj.p(:,1),...
                                             obj.p(:,2),obj.b,z,N,lon,lat);                   
end

%% Getting the J stuff if required (nonlocal tensor type)
if strcmp(type,'nonlocal')
   % Compute gradients of J from Nb, Nm and bathymetry B, and grid  
   [~,dJ] = Calc_dJ_Nyc_Unstruc(obj.t,obj.p,obj.b,Nm,omega,...
                                2,MinDepth,proj,[],4);                    
end

%% Calculate F_it from Nb, Nm, Jx, Jy, and H2_mesh or as reqd.
C_it = C_it / (4*pi); % Divide by 4pi here.
if isempty(N_filename)
    F_it = C_it * ones(size(obj.b));
    if strcmp(type(1:5),'local')
        H2_mesh = obj.bx.^2 + obj.by.^2;
    end
else
    Nb_t = Nb;
    if strcmp(crit,'Nmw')
        Nb_t = Nmw;
    end
    % I made C_it equivalent for directional, scalar and tensor (should be
    % about 0.20 for optimal)
    if strcmp(type(1:5),'local')
       F_it = C_it * sqrt((Nb_t.^2 - omega^2).*(Nm.^2 - omega^2))/omega;
    elseif strcmp(type,'nonlocal')
       F_it = C_it * Nb_t.*sqrt(1-f.^2/omega^2)./obj.b;
    end
    F_it = real(F_it); % in case becomes complex due to minus square root
    %
    H2_mesh = obj.bx.^2 + obj.by.^2;
    %% Compute criticality and normalise the friction
    alpha2  = (omega^2 - f.^2)./(Nb_t.^2 - omega^2);
    if ~strcmp(crit,'none')
        gamma2 = max(H2_mesh./alpha2,1);
        % Normalise F_it by criticality
        F_it = F_it./gamma2;
        % Make sure that if alpha2 < 0 that F_it is 
        % set equal to 0 since real part would be 0.
        F_it(alpha2 < 0) = 0;
    end
end

% Cut it off at the MinDepth
F_it(obj.b < MinDepth) = 0;

% Make zero for no slopes
F_it(obj.bx == 0 & obj.by == 0) = 0;
if strcmp(type,'local_sca')
    F_it(H2_mesh == 0) = 0;
elseif strcmp(type,'nonlocal')
    F_it(dJ(:,1) == 0 & dJ(:,1) == 0) = 0;
end
%
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
        if strcmp('internal_tide_friction',obj.f13.defval.Atr(NA).AttrName)
            broken = 1;
            % overwrite existing internal_tide
            break
        end
    end
    if ~broken 
        % add internal_tide to list
        obj.f13.nAttr = obj.f13.nAttr + 1;
        NA = obj.f13.nAttr;
    end
end

attrname = 'internal_tide_friction';  

% Default Values
obj.f13.defval.Atr(NA).AttrName = attrname;  
% We can just put in the options here
obj.f13.defval.Atr(NA).Unit = strcat(type,', C_it = ',num2str(C_it),...
                               ', ',num2str(crit),', D',num2str(MinDepth)); 
if isempty(N_filename)
    % We want to just give it the gradients of bathy (and J for the tensor)
    % multiplied by C_it and calculate the full values inside of ADCIRC
    % using N from say HYCOM
    if strcmp(type,'nonlocal')
        valpernode = 5;
    elseif strcmp(type,'local_dir')
        valpernode = 3;
    elseif strcmp(type,'local_sca')
        valpernode = 1;
    end
else
    % Provide the full coefficients
    if strcmp(type,'local_sca')
        valpernode = 1;
    else
        valpernode = 3;
    end
end   
obj.f13.defval.Atr(NA).ValuesPerNode = valpernode ;
obj.f13.defval.Atr(NA).Val = zeros(1,valpernode) ;

% User Values
obj.f13.userval.Atr(NA).AttrName = attrname;
numnodes = length(find(F_it > 0));
obj.f13.userval.Atr(NA).usernumnodes = numnodes ;
% Print out list of nodes for each
K = find(F_it > 0);
if isempty(N_filename)
    % We want to just give it the gradients of bathy (and J for the tensor)
    % multiplied by C_it and calculate the full values inside of ADCIRC
    % using N from say HYCOM
    if strcmp(type,'nonlocal')
        obj.f13.userval.Atr(NA).Val = [K F_it(K) obj.bx(K) ...
                                       obj.by(K) obj.dJ(K,:)]';
    elseif strcmp(type,'local_dir')
        obj.f13.userval.Atr(NA).Val = [K F_it(K) obj.bx(K) obj.by(K)]';
    elseif strcmp(type,'local_sca')
        obj.f13.userval.Atr(NA).Val = [K F_it(K) H2_mesh(K)]';
    end
else
    % Provide the full coefficients
    if strcmp(type,'local_sca')
        obj.f13.userval.Atr(NA).Val = [K F_it(K).*H2_mesh(K)]';
    elseif strcmp(type,'local_dir')
        obj.f13.userval.Atr(NA).Val = [K F_it(K).*(obj.bx(K).^2) ...
                     F_it(K).*(obj.by(K).^2) F_it(K).*(obj.bx(K).*obj.by(K))]';
    elseif strcmp(type,'nonlocal')
        obj.f13.userval.Atr(NA).Val = [K ...
                                    2*F_it(K).*(obj.bx(K).*obj.dJ(K,1)) ...
                                    2*F_it(K).*(obj.by(K).*obj.dJ(K,2)) ...
               F_it(K).*(obj.dJ(K,1).*obj.by(K) + obj.dJ(K,2).*obj.bx(K))];
        % Check if neg eigenvalue change to positive
        C = obj.f13.userval.Atr(NA).Val(:,2:end);
        for ii = 1:length(K)
            A = [C(ii,1) C(ii,3); C(ii,3) C(ii,2)];
            A = Eigen_neg_to_pos(A);
            obj.f13.userval.Atr(NA).Val(ii,2:end) = [A(1,1) A(2,2) A(2,1)];
        end
        obj.f13.userval.Atr(NA).Val = obj.f13.userval.Atr(NA).Val';
    end  
end 

if ~isempty(obj.f15)
    % Change attribute in obj.f15
    disp('Adding on itfric attribute name in fort.15 struct')
    obj.f15.nwp = obj.f15.nwp + 1;
    obj.f15.AttrName(obj.f15.nwp).name = attrname;
end

%EOF
end