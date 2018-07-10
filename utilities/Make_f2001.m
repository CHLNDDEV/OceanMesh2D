function obj = Make_f2001(obj,File_Prefix,File_Suffixes,bathyfile,ETIMINC)
% obj = Make_f2001(obj,File_Prefix,File_Suffixes,bathyfile,ETIMINC)
% Input a msh class object get the values of the elevations and velocities
% to put in the sponge for times between ts and te based on the 
% File_Suffixes of the  input files (location is File_Prefix). 
% Bathyfile is the file of bathymetric values used to calculate fluxes from
% the velocities before dividing by the msh class depths to get back to
% velocities
%
%  Author:      William Pringle                                 
%  Created:     May 15 2018                                      
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

BZ = zeros(length(sp_x),length(File_Suffixes));
BFX = zeros(length(sp_x),length(File_Suffixes));
BFY = zeros(length(sp_x),length(File_Suffixes));
for t = 1:length(File_Suffixes)
    filename = strcat(File_Prefix,File_Suffixes(t),".nc");
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

        % Delete by range search
        Krs = rangesearch([lon_x,lat_y],sp_points,0.25); % 0.25 deg radius
        LKd = 0; Kd = [];
        for i = 1:length(sp_points)
           K = Krs{i};
           J = LKd + 1:LKd + length(K);
           Kd(J) = K; 
           LKd = length(Kd);
        end
        Kd = unique(Kd);

        % The new lon and lat vectors of data
        lon_x = lon_x(Kd); lat_y = lat_y(Kd);

        % Do the projection
        [x,y] = m_ll2xy(lon_x,lat_y);       
    end
    
    % surface elevations
    zeta = ncread(filename,'surf_el');
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
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),qx(~isnan(Z)),'natural');
    BFX(:,t) = F(sp_x,sp_y)./sp_b;  
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),qy(~isnan(Z)),'natural');
    BFY(:,t) = F(sp_x,sp_y)./sp_b;  
end

%% Make into f2001 struct
filename1 = strcat(File_Prefix,File_Suffixes(1),".nc");
filename2 = strcat(File_Prefix,File_Suffixes(end),".nc");
[~,n1] = fileparts(char(filename1)); [~,n2] = fileparts(char(filename2));
obj.f2001.DataTitle = [n1 ' -> ' n2];
obj.f2001.ETIMINC = ETIMINC;
obj.f2001.NumOfNodes = size(BZ,1);
obj.f2001.NumOfTimes = size(BZ,2);
obj.f2001.Val = [BZ(:) BFX(:) BFY(:)];
%EOF
end


