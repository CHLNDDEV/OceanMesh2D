function obj = GridData(geodata,obj,varargin)
% obj = GridData(geodata,obj,'K',K,'type',type);
% GridData: Uses the cell-averaged approach to interpolate the nodes
%           on the unstructured grid to the DEM.          
%       Input : geodata - either a geodata class or a filename of the netcdf DEM
%                   obj - A msh class object containing:
%                     p - A vector of the locations of all the vertices
%                         in the mesh, length in number of nodes, nn x 2
%                     t - Matrix size ne x 3 of the element table of the mesh
%
%          K (optional) - vector of relevant nodes to search. This can
%                         significantly speed up the calculation by either 
%                         only interpolating part of the mesh or 
%                         intepolating whole mesh but in segments. Function
%                         automatically removes uncessary portions of DEM
%                         Example to call: 
%                         K = find( obj.p(:,1) >= lon_min & ...
%                                   obj.p(:,2) <= lon_max);
%                         obj = GridData(fname,obj,'K',K);
%       type (optional) - type is either 'depth', 'slope' or 'all'
%                         'all' is the default (both slope and depth). 
%                         'slope' to gets the gradients of DEM
%
%     interp (optional) - interp is either the normal griddedInterpolant
%                         options in MATLAB or is 'CA' (default). Note: 'CA'
%                         applies linear griddedInterpolant when the DEM
%                         and grid sizes are similar. 
%
%          N (optional) - enlarge cell-averaging stencil by factor N (only 
%                         relevant for CA interpolation method). 
%                         default value N=1. 
%
%        nan (optional) - 'fill' to fill in any NaNs appearing in bathy
%
%   mindepth (optional) - ensure the minimum depth is bounded in the 
%                         interpolated region 
%
%   maxdepth (optional) - ensure the maximum depth is bounded in the 
%                         interpolated region 
%
%   ignoreOL (optional) - NaN overland data for more accurate seabed interpolation
%
%      Output : obj     - A mesh class object with the following updated:
%               b       - The depth or vertical difference from the datum
%                         in the mesh if type is 'depth' or 'all'
%               bx      - The x direction gradient if type is 'slope' or
%                         'all'
%               by      - The y direction gradient if type is 'slope' or
%                         'all'
%
%	Author: William Pringle, and Keith Roberts CHL, Notre Dame University
%	Created: 2018-03-12
%
%   Edits by Keith Roberts, 2018-02-22 to avoid FillValue contaimation and
%   to only interpolate data inside extent of DEM. 
%
%   Edits by Keith Roberts, 2018-10-18 to bound depth in interpolated
%   regions
%
%   Edits by Keith Roberts, 2019-4-4 to interpolate the floodplain

%% Test optional arguments
% default
nn = length(obj.p);
K  = (1:nn)';        % K is all of the grid
type = 'all';
interp = 'CA';
NaNs = 'ignore';
N    = 1 ; 
mindepth = -inf ; 
maxdepth = +inf ; 
if ~isempty(varargin)
    varargin=varargin{1} ; 
    names = {'K','type','interp','nan','N','mindepth','maxdepth','ignoreOL'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                K = varargin{ind*2}; 
                nn = length(K);
            elseif ii == 2
                type = varargin{ind*2};
            elseif ii == 3
                interp = varargin{ind*2};
            elseif ii == 4
                NaNs = varargin{ind*2};
            elseif ii == 5 
                N = varargin{ind*2} ; 
            elseif ii ==6 
              mindepth = varargin{ind*2} ; 
            elseif ii ==7 
              maxdepth = varargin{ind*2}  ; 
            elseif ii ==8
              ignoreOL = varargin{ind*2} ; 
            end
        end    
    end
end

if ~exist('ignoreOL','var')
  ignoreOL = 0 ; 
end
if ignoreOL 
  disp('NaNing overland data before interpolating') 
end
if N > 1 
   disp(['Enlarging CA stencil by factor ',num2str(N)]) ;  
end
if mindepth > -inf 
   disp(['Bounding minimum depth to ',num2str(mindepth), ' meters.']) ;  
end
if maxdepth < inf 
   disp(['Bounding maximum depth to ',num2str(maxdepth), ' meters.']) ;  
end

if strcmp(type,'slope')
    warning(['You must have the bathymetry on the grid before '...
            'calculating the slopes'])
end

if strcmp(NaNs,'fill')
   disp('Fill in NaNs using nearest neighbour interpolation.')
   warning(['Note that will try and put bathy everywhere on mesh even ' ...
           'outside of gdat/dem extents unless K logical is set.'])
end

%% Let's read the LON LAT of DEM if not already geodata
flipUD = 0;
if ~isa(geodata,'geodata')
    try
        DEM_XA = double(ncread(geodata,'lon'));
        DEM_YA = double(ncread(geodata,'lat'));
    catch
        DEM_XA = double(ncread(geodata,'x'));
        DEM_YA = double(ncread(geodata,'y'));
    end
    DELTA_X = mean(diff(DEM_XA));
    DELTA_Y = mean(diff(DEM_YA));
    if DELTA_Y < 0
       flipUD = 1;
       DEM_YA = flipud(DEM_YA);
       DELTA_Y = mean(diff(DEM_YA));
    end
else
    DEM_XA = geodata.Fb.GridVectors{1};
    DEM_YA = geodata.Fb.GridVectors{2};
    DEM_Z  = -geodata.Fb.Values;
    DELTA_X = mean(diff(DEM_XA));
    DELTA_Y = mean(diff(DEM_YA));
    [DEM_X,DEM_Y] = ndgrid(DEM_XA,DEM_YA);
end

if max(DEM_XA) > 180
   lon_change = obj.p(:,1) < 0; 
else
   lon_change = false(length(obj.p),1);
end
obj.p(lon_change,1) = obj.p(lon_change,1) + 360;

% kjr edit 20180320
if length(K) == length(obj.p)
    % msh class may not have populated b, if this case populate with NaN to
    % detect issues
    if isempty(obj.b)
        obj.b = obj.p(:,1)*NaN; 
    end
    if ~strcmp(NaNs,'fill')
        outside = obj.p(:,1) < min(DEM_XA)-DELTA_X | ...
                  obj.p(:,1) > max(DEM_XA)+DELTA_X | ...
                  obj.p(:,2) < min(DEM_YA)-DELTA_Y | ...
                  obj.p(:,2) > max(DEM_YA)+DELTA_Y;
        K(outside) = [];
    end
end

% If the length of DEM_X and DEM_Y is too large then let's break it up by
% using K 
% test 
lon_min = min(obj.p(K,1));
lon_max = max(obj.p(K,1));
lat_min = min(obj.p(K,2));
lat_max = max(obj.p(K,2));

I = find(DEM_XA >= lon_min & DEM_XA <= lon_max);
J = find(DEM_YA >= lat_min & DEM_YA <= lat_max);
% single or integer array
L = length(I)*length(J)*4e-9;
% larger than 2 GB
if L > 2
    % divide by latitude
    times = ceil(L/2); lat_l = lat_min; dl = (lat_max - lat_min)/times;
    for t = 1:times
        lat_r = lat_l + dl;
        if t == times
            lat_r = lat_max;
        end
        KK{t} = K(obj.p(K,2) >= lat_l & obj.p(K,2) <= lat_r);
        lat_l = lat_r;
    end
else
    times = 1;
    KK{1} = K;
end
% This deletes any elements straddling the -180/180 boundary 
xt = [obj.p(obj.t(:,1),1) obj.p(obj.t(:,2),1) ...
      obj.p(obj.t(:,3),1) obj.p(obj.t(:,1),1)];
dxt = diff(xt,[],2);
tt = obj.t;
tt(abs(dxt(:,1)) > 180 | abs(dxt(:,2)) > 180 |  abs(dxt(:,2)) > 180,:) = [];

% Do this once;
vtoe_o = VertToEle(tt); %find connecting elements to each node
pmid = squeeze(mean(reshape(obj.p(tt,:),[],3,2),2)); % get mid points of elements
pmid(end+1,:) = NaN;

K_o = K; % save the original K
disp(['Looping over the DEM ' num2str(times) ' time(s)'])
for t = 1:times
K = KK{t};
%% Get the grid size for each node
% Get the subset 
vtoe = vtoe_o(:,K);
% make sure when vtoe is zero we just return a NaN
vtoe(vtoe == 0) = length(tt) + 1;
% the connecting element centers
pcon = reshape(pmid(vtoe,:),size(vtoe,1),[],2);

% delta is max and min bounds of the connecting element centers (or if only
% one connecting element then do difference with the vertex itself

% kjr 10/17/2018 enlarge the cell-averaging stencil by a factor N 
pcon_max = max(squeeze(max(pcon,[],1)),2*obj.p(K,:)-squeeze(min(pcon,[],1)));
pcon_min = min(squeeze(min(pcon,[],1)),2*obj.p(K,:)-squeeze(max(pcon,[],1)));
delta = pcon_max - pcon_min;

pcon_max = N*pcon_max+(1-N)*obj.p(K,:) ;
pcon_min = N*pcon_min+(1-N)*obj.p(K,:) ;
                

%% Now read the bathy (only what we need)
if ~isa(geodata,'geodata')
    BufferL = max(delta); 
    lon_min = min(obj.p(K,1)) - BufferL(1);
    lon_max = max(obj.p(K,1)) + BufferL(1);
    lat_min = min(obj.p(K,2)) - BufferL(2);
    lat_max = max(obj.p(K,2)) + BufferL(2);
    I = find(DEM_XA > lon_min & DEM_XA < lon_max);
    J = find(DEM_YA > lat_min & DEM_YA < lat_max);
    if exist('DEM_X','var')
        clear DEM_X DEM_Y DEM_Z
    end
    [DEM_X,DEM_Y] = ndgrid(DEM_XA(I),DEM_YA(J));
    
    finfo = ncinfo(geodata);
    for ii = 1:length(finfo.Variables)
        if length(finfo.Variables(ii).Size) == 2
            zvarname = finfo.Variables(ii).Name;
            break
        end
    end
    DEM_Z = single(ncread(geodata,zvarname,...
                   [I(1) J(1)],[length(I) length(J)]));
    if flipUD
       DEM_Z = fliplr(DEM_Z);
    end
               
    % make into depths (ADCIRC compliant) 
    DEM_Z = -DEM_Z;
end

% bound all depths below mindepth 
DEM_Z(DEM_Z < mindepth) = mindepth ; 

% bound all depths above maxdepth 
DEM_Z(DEM_Z > maxdepth) = maxdepth ; 

%% Make the new bx, by, b and calculate gradients if necessary
if strcmp(type,'slope') || strcmp(type,'all')
    bx = NaN(length(K),1); 
    by = NaN(length(K),1); 
    % Get the dx and dy of the dem in meters
    DELTA_X1 = DELTA_X*111e3;
    DELTA_X1 = DELTA_X1*cosd(DEM_Y(1,:)); % for gradient function   
    DELTA_Y1 = DELTA_Y*111e3;
    %m_idist(mean(obj.p(K,1)),mean(obj.p(K,2)),...
    %        mean(obj.p(K,1)),mean(obj.p(K,2))+DELTA_Y)
    %m_idist(mean(obj.p(K,1)),mean(obj.p(K,2)),...
    %        mean(obj.p(K,1))+DELTA_X,mean(obj.p(K,2)))
    if exist('DEM_ZY','var')
        clear DEM_ZY DEM_ZX
    end
    % Calculate the gradient of the distance function.
    [DEM_ZY,DEM_ZX] = EarthGradient(DEM_Z,DELTA_Y1,DELTA_X1);
    %[DEM_ZY,DEM_ZX] = gradient(DEM_Z,DELTA_Y1,mean(DELTA_X1));
    % New method of averaging the absolute values 
    DEM_ZY = abs(DEM_ZY); DEM_ZX = abs(DEM_ZX);
end
if strcmp(type,'depth') || strcmp(type,'all')
    b = NaN(length(K),1); 
end

%% Calculate N - number of DEM grids to average - for each node
tic
if strcmp(interp,'CA')
    [~,IDXX,IDXY] = FindLinearIdx(obj.p(K,1),obj.p(K,2),DEM_X,DEM_Y);
    [~,IDXR,IDXT] = FindLinearIdx(pcon_max(:,1),pcon_max(:,2),DEM_X,DEM_Y);
    [~,IDXL,IDXB] = FindLinearIdx(pcon_min(:,1),pcon_min(:,2),DEM_X,DEM_Y);
    % make sure inside box
    % y
    high = DEM_Y(IDXR(1),IDXT)' > pcon_max(:,2);
    IDXT(high) = max(1,IDXT(high) - 1);
    low = DEM_Y(IDXL(1),IDXB)' < pcon_min(:,2);
    IDXB(low) = min(size(DEM_Y,2),IDXB(low) + 1);
    % x
    high = DEM_X(IDXR,IDXT(1)) > pcon_max(:,1);
    IDXR(high) = max(1,IDXR(high) - 1);
    low = DEM_X(IDXL,IDXB(1)) < pcon_min(:,1);
    IDXL(low) = min(size(DEM_X,1),IDXL(low) + 1);
    % Make sure no negative or positive Nx, Ny
    I = IDXL > IDXR;
    IDXL(I) = IDXX(I); IDXR(I) = IDXX(I);
    I = IDXB > IDXT;
    IDXB(I) = IDXY(I); IDXT(I) = IDXY(I);
    % The span of x and y
    Nx = IDXR - IDXL;
    Ny = IDXT - IDXB;
    
    disp(['Performing cell-averaging with Ny(max) ' num2str(max(Ny)) ...
          ' and Nx(max)' num2str(max(Nx))])
      
    % Average for the depths    
    if strcmp(type,'depth') || strcmp(type,'all')
        for ii = 1:length(K)
            %%% MTC special option goes here
            if(ignoreOL)
              pts = reshape(DEM_Z(IDXL(ii):IDXR(ii),...
                  IDXB(ii):IDXT(ii)),[],1);
              pts(pts < 0) = NaN ;
            else
              pts = reshape(DEM_Z(IDXL(ii):IDXR(ii),...
                  IDXB(ii):IDXT(ii)),[],1);
           end
           b(ii) = mean(pts,'omitnan');            
        end
        % Try and fill in the NaNs
        if strcmp(NaNs,'fill')
            if ~isempty(find(isnan(b),1))
                localcoord = obj.p(K,:);
                KI = knnsearch(localcoord(~isnan(b),:),localcoord(isnan(b),:));
                bb = b(~isnan(b),:);
                b(isnan(b)) = bb(KI); clear bb localcoord
            end
        end
    end
        
    % Average for the slopes
    if strcmp(type,'slope') || strcmp(type,'all')
        for ii = 1:length(K)
            bx(ii) = mean(reshape(DEM_ZX(IDXL(ii):IDXR(ii),...
                                         IDXB(ii):IDXT(ii)),[],1),'omitnan');
            by(ii) = mean(reshape(DEM_ZY(IDXL(ii):IDXR(ii),...
                                         IDXB(ii):IDXT(ii)),[],1),'omitnan');
        end
        if strcmp(NaNs,'fill')
            % Try and fill in the NaNs
            if ~isempty(find(isnan(bx),1))
                localcoord = obj.p(K,:);
                KI = knnsearch(localcoord(~isnan(bx),:),localcoord(isnan(bx),:));
                bb = bx(~isnan(bx),:);
                bx(isnan(bx)) = bb(KI); clear bb localcoord
            end
            if ~isempty(find(isnan(by),1))
                localcoord = obj.p(K,:);
                KI = knnsearch(localcoord(~isnan(by),:),localcoord(isnan(by),:));
                bb = by(~isnan(by),:);
                by(isnan(by)) = bb(KI); clear bb localcoord
            end
        end
    end
else
    %% Get gridded interpolant for non-CA interp
    disp(['Using griddedInterpolant with ' interp 'interp option'])
    if strcmp(type,'slope') || strcmp(type,'all')
        Fx = griddedInterpolant(DEM_X,DEM_Y,DEM_ZX,interp);
        Fy = griddedInterpolant(DEM_X,DEM_Y,DEM_ZY,interp);
        bx = Fx(obj.p(K,1),obj.p(K,2));
        by = Fy(obj.p(K,1),obj.p(K,2));
    end
    if strcmp(type,'depth') || strcmp(type,'all')
        F = griddedInterpolant(DEM_X,DEM_Y,DEM_Z,interp,'none');
        b = F(obj.p(K,1),obj.p(K,2));
    end
end
toc
%% Put in the global msh class
if strcmp(type,'depth') || strcmp(type,'all')
    if isempty(obj.b)
        obj.b = zeros(length(obj.p),1);
    end
    % Only overwrite non-nan values
    obj.b(K(~isnan(b))) = b(~isnan(b));
end
if strcmp(type,'slope') || strcmp(type,'all')
    if isempty(obj.bx)
        obj.bx = zeros(length(obj.p),1);
    end 
    obj.bx(K(~isnan(bx))) = bx(~isnan(bx));
    % Put in the global array
    if isempty(obj.by)
        obj.by = zeros(length(obj.p),1);
    end
    obj.by(K(~isnan(by))) = by(~isnan(by));
end
end
% New method of averaging asbolute values of slope before multiplying
% by the sign of the slope on unstructured grid
if strcmp(type,'slope') || strcmp(type,'all')
    [Hx,Hy] = Unstruc_Bath_Slope( tt,obj.p(:,1),obj.p(:,2),obj.b);
    obj.bx(K_o) = sign(Hx(K_o)).*obj.bx(K_o); 
    obj.by(K_o) = sign(Hy(K_o)).*obj.by(K_o);
end
obj.p(lon_change,1) = obj.p(lon_change,1) - 360;
%EOF
end
