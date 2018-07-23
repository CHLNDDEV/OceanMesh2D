function obj = Make_f2001(obj,filename,bathyfile)
% obj = Make_f2001(obj,filename,bathyfile)
% Input a msh class object get the values of the elevations and velocities
% to put in the sponge for times in the filename.
% Bathyfile is the file of bathymetric values used to calculate fluxes from
% the velocities before dividing by the msh class depths to get back to
% velocities
%
%  Author:      William Pringle                                 
%  Created:     May 15 2018, 
%               July 19 2018, Updated to do it just for one file                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = find(contains({obj.f13.defval.Atr(:).AttrName},'sponge'));
if isempty(ii)
    error('No Sponge infomation to use')
end

%% Get sponge info
nodes = obj.f13.userval.Atr(ii).Val(1,:)';
sp_points = obj.p(nodes,:);
sp_b = max(10,obj.b(nodes)); % make sure not less than 10 m

%% Do some projection for sp_points
proj = 'Mercator';
m_proj(proj,'lon',[ min(sp_points(:,1))-0.25 max(sp_points(:,1))+0.25 ],...
            'lat',[ min(sp_points(:,2))-0.25 max(sp_points(:,2))+0.25]) 

% Do projection
[sp_x,sp_y] = m_ll2xy(sp_points(:,1),sp_points(:,2));   

time = ncread(filename,'time');
DT = median(diff(time))*3600; %[s]
L = length(time);
BZ = zeros(length(sp_x),L);
BFX = zeros(length(sp_x),L);
BFY = zeros(length(sp_x),L);
disp(['Interpolating ' num2str(L) ' snaps'])
for t = 1:L
    if mod(t,100) == 0 
        disp(['Snap #' num2str(t)])
    end
    %% Read the grid and density data 
    if t == 1
        % read lon, lat and standard depths
        lon = ncread(filename,'lon');
        lat = ncread(filename,'lat');
        z = ncread(filename,'depth');
        lon(lon > 180) = lon(lon > 180) - 360;
        
        % get local bathymetry
        b = ncread(bathyfile,'depth');
        lont = ncread(bathyfile,'Longitude');
        lont(lont > 180) = lont(lont > 180) - 360;
        latt = ncread(bathyfile,'Latitude');
        [~,IA] = intersect(lont,lon);
        [~,IB] = intersect(latt,lat);
        b = b(IA,IB);
        
        [lon,lat] = ndgrid(lon,lat);
        lon_x = reshape(lon,[],1);
        lat_y = reshape(lat,[],1);
        % Delete uncessecary portions
        % First delete by square
        I = find(lon_x < min(sp_points(:,1))-0.25 | ...
                 lon_x > max(sp_points(:,1))+0.25 | ...
                 lat_y < min(sp_points(:,2))-0.25 | ...
                 lat_y > max(sp_points(:,2))+0.25);
        lon_x(I) = []; lat_y(I) = []; 

        % Delete by our knnsearch
        Kd = ourKNNsearch([lon_x,lat_y]',sp_points',20);
        Kd = unique(Kd);
        
        % The new lon and lat vectors of data
        lon_x = lon_x(Kd); lat_y = lat_y(Kd);

        % Do the projection
        [x,y] = m_ll2xy(lon_x,lat_y);       
    end
    
    % surface elevations
    zeta = ncread(filename,'surf_el',[1 1 t],[size(lon) 1]);
    Z = reshape(zeta,[],1);
    Z(I) = []; Z = Z(Kd);
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),Z(~isnan(Z)),'natural');
    BZ(:,t) = F(sp_x,sp_y);  
    
    % fluxes
    u = ncread(filename,'water_u');
    u(isnan(u)) = 0;
    v = ncread(filename,'water_v');
    v(isnan(v)) = 0;
    % this loop makes sure we only take the integral over the actual depth
    % of the source mesh
    qx = 0*b;  qy = 0*b; 
    for ii = 1:size(b,1)
        for jj = 1:size(b,2)
            zn = z;
            zn(z > b(ii,jj)) = b(ii,jj);       
            qx(ii,jj) = trapz(zn,squeeze(u(ii,jj,:)));
            qy(ii,jj) = trapz(zn,squeeze(v(ii,jj,:)));
        end
    end
    qx = reshape(qx,[],1);
    qx(I) = []; qx = qx(Kd);
    qy = reshape(qy,[],1);
    qy(I) = []; qy = qy(Kd);
    % do interpolation and divide by our own depths to get velocities
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),...
                             qx(~isnan(Z)),'natural');
    BFX(:,t) = F(sp_x,sp_y)./sp_b;  
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),...
                             qy(~isnan(Z)),'natural');
    BFY(:,t) = F(sp_x,sp_y)./sp_b;  
end

%% Make into f2001 struct
ts = datenum('2000-01-01') + time(1)/24;
te = datenum('2000-01-01') + time(end)/24;
obj.f2001.DataTitle = [datestr(ts) ' -> ' datestr(te)];
obj.f2001.ETIMINC = DT;
obj.f2001.NumOfNodes = size(BZ,1);
obj.f2001.NumOfTimes = size(BZ,2);
obj.f2001.Val = [BZ(:) BFX(:) BFY(:)];
%EOF
end


