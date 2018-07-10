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
%       type (optional) - type is either 'depth', 'slope' or 'all'. 
%                         'all' is the default (both slope and depth). 
%                         'slope' to gets the gradients of DEM
%
%     interp (optional) - interp is either the normal griddedInterpolant
%                         options in MATLAB or is 'CA' (default). Note: 'CA'
%                         applies linear griddedInterpolant when the DEM
%                         and grid sizes are similar. 
%
%      Output : obj     - A mesh class object with the following updated:
%               b       - The depth or vertical difference from the datum
%                         in the mesh if type is 'depth' or 'all'
%               bx      - The x direction gradient if type is 'slope' or
%                         'all'
%               by      - The y direction gradient if type is 'slope' or
%                         'all'
%
%	Author: William Pringle, CHL, Notre Dame University
%	Created: 2018-03-12
%   Edits by Keith Roberts, 2018-02-22 to avoid FillValue contaimation and
%   to only interpolate data inside extent of DEM. 

%% Test optional arguments
% default
nn = length(obj.p);
K  = (1:nn)';        % K is all of the grid
type = 'all';
interp = 'CA';
if ~isempty(varargin)
    varargin=varargin{1} ; 
    names = {'K','type','interp'};
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
            end
        end    
    end
end

if strcmp(type,'slope')
    warning(['You must have the bathymetry on the grid before '...
            'calculating the slopes'])
end

%% Let's read the LON LAT of DEM if not already geodata
if ~isa(geodata,'geodata')
    DEM_XA = double(ncread(geodata,'lon'));
    DEM_YA = double(ncread(geodata,'lat'));
    DELTA_X = mean(diff(DEM_XA));
    DELTA_Y = mean(diff(DEM_YA));
else
    DEM_XA = geodata.Fb.GridVectors{1};
    DEM_YA = geodata.Fb.GridVectors{2};
    DEM_Z  = -geodata.Fb.Values;
    DELTA_X = mean(diff(DEM_XA));
    DELTA_Y = mean(diff(DEM_YA));
    [DEM_X,DEM_Y] = ndgrid(DEM_XA,DEM_YA);
end

% kjr edit 20180320
if length(K) == length(obj.p)
    % msh class may not have populated b, if this case populate with NaN to
    % detect issues
    if isempty(obj.b)
        obj.b = obj.p(:,1)*NaN; 
    end
    outside = obj.p(:,1) <= min(DEM_XA(:)) | obj.p(:,1) >= max(DEM_XA(:)) | ...
              obj.p(:,2) <= min(DEM_YA(:)) | obj.p(:,2) >= max(DEM_YA(:));
    K(outside) = [];
end

% If the length of DEM_X and DEM_Y is too large then let's break it up by
% using the K
% test 
lon_min = min(obj.p(K,1));
lon_max = max(obj.p(K,1));
lat_min = min(obj.p(K,2));
lat_max = max(obj.p(K,2));

I = find(DEM_XA > lon_min & DEM_XA < lon_max);
J = find(DEM_YA > lat_min & DEM_YA < lat_max);
% single or integer array
L = length(I)*length(J)*4e-9;
% larger than 2 GB
if L > 2
    times = ceil(L/2); lon_l = lon_min; dl = (lon_max - lon_min)/times;
    for t = 1:times
        lon_r = lon_l + dl;
        if t == times
            lon_r = lon_max;
        end
        KK{t} = K(obj.p(K,1) >= lon_l & obj.p(K,1) <= lon_r);
        lon_l = lon_r;
    end
else
    times = 1;
    KK{1} = K;
end
% Do this once;
vtoe_o = VertToEle(obj.t); %find connecting elements to each node
pmid = squeeze(mean(reshape(obj.p(obj.t,:),[],3,2),2)); % get mid points of elements
pmid(end+1,:) = NaN;

K_o = K; % save the original K
disp(['Looping over the DEM ' num2str(times) ' time(s)'])
for t = 1:times
K = KK{t};
%% Get the grid size for each node
% Get the subset 
vtoe = vtoe_o(:,K);
% make sure when vtoe is zero we just return a NaN
vtoe(vtoe == 0) = length(obj.t) + 1;
% the connecting element centers
pcon = reshape(pmid(vtoe,:),size(vtoe,1),[],2);
% delta is max and min bounds of the connecting element centers (or if only
% one connecting element then do difference with the vertex itself
pcon_max = max(squeeze(max(pcon,[],1)),2*obj.p(K,:)-squeeze(min(pcon,[],1)));
pcon_min = min(squeeze(min(pcon,[],1)),2*obj.p(K,:)-squeeze(max(pcon,[],1)));
delta = pcon_max - pcon_min;

%% Now read the bathy (only what we need)
if ~isa(geodata,'geodata')
    BufferL = max(delta); 
    lon_min = min(obj.p(K,1)) - BufferL(1);
    lon_max = max(obj.p(K,1)) + BufferL(1);
    lat_min = min(obj.p(K,2)) - BufferL(2);
    lat_max = max(obj.p(K,2)) + BufferL(2);
    I = find(DEM_XA > lon_min & DEM_XA < lon_max);
    J = find(DEM_YA > lat_min & DEM_YA < lat_max);
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
    % make into depths
    DEM_Z = -DEM_Z;
end

%% Make the new bx, by, b and calculate gradients if necessary
if strcmp(type,'slope') || strcmp(type,'all')
    bx = NaN(length(K),1); 
    by = NaN(length(K),1); 
    % Get the dx and dy of the dem in meters
    DELTA_X1 = m_idist(mean(obj.p(K,1)),mean(obj.p(K,2)),...
                      mean(obj.p(K,1))+DELTA_X,mean(obj.p(K,2)));
    DELTA_Y1 = m_idist(mean(obj.p(K,1)),mean(obj.p(K,2)),...
                      mean(obj.p(K,1)),mean(obj.p(K,2))+DELTA_Y);             
    [DEM_ZY,DEM_ZX] = gradient(DEM_Z,DELTA_Y1,DELTA_X1);
    % New method of averaging the asbolute values 
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
    IDXT(high) = IDXT(high) - 1;
    low = DEM_Y(IDXL(1),IDXB)' < pcon_min(:,2);
    IDXB(low) = IDXB(low) + 1;
    % x
    high = DEM_X(IDXR,IDXT(1)) > pcon_max(:,1);
    IDXR(high) = IDXR(high) - 1;
    low = DEM_X(IDXL,IDXB(1)) < pcon_min(:,1);
    IDXL(low) = IDXL(low) + 1;
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
            b(ii) = mean(reshape(DEM_Z(IDXL(ii):IDXR(ii),IDXB(ii):IDXT(ii)),[],1),'omitnan');
        end
    end
    % Average for the slopes
    if strcmp(type,'slope') || strcmp(type,'all')
        for ii = 1:length(K)
            bx(ii) = mean(reshape(DEM_ZX(IDXL(ii):IDXR(ii),IDXB(ii):IDXT(ii)),[],1),'omitnan');
            by(ii) = mean(reshape(DEM_ZY(IDXL(ii):IDXR(ii),IDXB(ii):IDXT(ii)),[],1),'omitnan');
        end
    end
%     b1 = b; bx1 = bx; by1 = by;
%     tic
%     ny = size(DEM_X,1); nx = size(DEM_X,2);            
%     Nx_un = unique(Nx);
%     for Nxx = Nx_un' 
%         Ny_un = unique(Ny(Nx == Nxx));
%         for Nyy = Ny_un' 
%             [ipos,jpos] = ind2sub(size(DEM_X),IDX(Nx == Nxx & Ny == Nyy));
%             npos = zeros((Nxx+1)*(Nyy+1),length(ipos));
%             nn = 0;
%             for i = -ceil(Nyy/2):ceil(Nyy/2)
%                 for j = -ceil(Nxx/2):ceil(Nxx/2)
%                     nn = nn + 1;
%                     npos(nn,:) = (jpos+j)*ny + max(1,min(ipos+i,nx));
%                 end
%             end     
%             % make sure the values are not outside the DEM
%             npos(npos <= 0 | npos > numel(DEM_X)) = length(DEM_Z);
%             if strcmp(type,'depth') || strcmp(type,'all')
%                 b(Nx == Nxx & Ny == Nyy) = mean(DEM_Z(npos));
%             end
%             if strcmp(type,'slope') || strcmp(type,'all')
%                 bx(Nx == Nxx & Ny == Nyy) = mean(DEM_ZX(npos));
%                 by(Nx == Nxx & Ny == Nyy) = mean(DEM_ZY(npos));
%             end
%         end
%     end
%     toc
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
        F = griddedInterpolant(DEM_X,DEM_Y,DEM_Z,interp);
        b = F(obj.p(K,1),obj.p(K,2));
    end
end
toc
%% Put in the global msh class
if strcmp(type,'depth') || strcmp(type,'all')
    if isempty(obj.b)
        obj.b = zeros(length(obj.p),1);
    end
    obj.b(K) = b;
end
if strcmp(type,'slope') || strcmp(type,'all')
    if isempty(obj.bx)
        obj.bx = zeros(length(obj.p),1);
    end 
    obj.bx(K) = bx;
    % Put in the global array
    if isempty(obj.by)
        obj.by = zeros(length(obj.p),1);
    end
    obj.by(K) = by;
end
end
% New method of averaging asbolute values of slope before multiplying
% by the sign of the slope on unstructured grid
if strcmp(type,'slope') || strcmp(type,'all')
    [Hx,Hy] = Unstruc_Bath_Slope( obj.t,obj.p(:,1),obj.p(:,2),obj.b);
    obj.bx(K_o) = sign(Hx(K_o)).*obj.bx(K_o); 
    obj.by(K_o) = sign(Hy(K_o)).*obj.by(K_o);
end
%EOF
end