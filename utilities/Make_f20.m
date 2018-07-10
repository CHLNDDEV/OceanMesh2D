function obj = Make_f20(obj,File_Prefix,File_Suffixes,bathyfile,ETIMINC)
% obj = Make_f20(obj,File_Prefix,File_Suffixes,bathyfile,ETIMINC)
% Input a msh class object get the values of the elevations and normal
% fluxes for times between ts and te based on the File_Suffixes of the 
% input files (location is File_Prefix). Bathyfile is the file of 
% bathymetric values used to calculate fluxes from the velocities
%
%  Author:      William Pringle                                 
%  Created:     May 7 2018                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(obj.op)
    error('No boundary information to interpolate to')
end

%% Do some projection for obj.p
proj = 'Mercator';
m_proj(proj,'lon',[ min(obj.p(:,1))-0.25 max(obj.p(:,1))+0.25 ],...
            'lat',[ min(obj.p(:,2))-0.25 max(obj.p(:,2))+0.25]) 

%% Get boundary info
opedat = obj.op;
b_lon = zeros(opedat.neta,1);
b_lat = zeros(opedat.neta,1);
node_num = zeros(opedat.neta,1);
ns = 1; ne = 0;
for n = 1:opedat.nope
    ne = ne + opedat.nvdll(n);
    node_num(ns:ne) = opedat.nbdv(1:opedat.nvdll(n),n);
    b_lon(ns:ne)    = obj.p(node_num(ns:ne),1);
    b_lat(ns:ne)    = obj.p(node_num(ns:ne),2);
    ns = ne + 1;
end
% Do projection
[b_x,b_y] = m_ll2xy(b_lon,b_lat);   

BZ = zeros(length(b_x),length(File_Suffixes));
BF = zeros(length(b_x),length(File_Suffixes));
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
        I = find(lon_x < min(obj.p(:,1))-0.25 | lon_x > max(obj.p(:,1))+0.25 | ...
                 lat_y < min(obj.p(:,2))-0.25 | lat_y > max(obj.p(:,2))+0.25);
        lon_x(I) = []; lat_y(I) = []; 

        % Delete by range search
        Krs = rangesearch([lon_x,lat_y],[b_lon, b_lat],0.25); % 0.25 deg radius
        LKd = 0; Kd = [];
        for i = 1:length(b_lon)
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
    BZ(:,t) = F(b_x,b_y);  
    
    % fluxes
    u = ncread(filename,'water_u');
    u(isnan(u)) = 0;
    % this loop makes sure we only take the integral over the actual depth
    % of the source mesh
    qx = 0*b;
    for ii = 1:size(b,1)
        for jj = 1:size(b,2)
            zn = z;
            zn(z > b(ii,jj)) = b(ii,jj);       
            qx(ii,jj) = trapz(zn,squeeze(u(ii,jj,:)));
        end
    end
    % qx = trapz(z,u,3); %this line will do the integration based on the 
    % standard depth contours
    qx = reshape(qx,[],1);
    qx(I) = []; qx = qx(Kd);
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),qx(~isnan(Z)),'natural');
    BF(:,t) = F(b_x,b_y);  
end

%% Make into f20 struct
filename1 = strcat(File_Prefix,File_Suffixes(1),".nc");
filename2 = strcat(File_Prefix,File_Suffixes(end),".nc");
[~,n1] = fileparts(char(filename1)); [~,n2] = fileparts(char(filename2));
obj.f20.DataTitle = [n1 ' -> ' n2];
obj.f20.ETIMINC = ETIMINC;
obj.f20.NumOfNodes = size(BZ,1);
obj.f20.NumOfTimes = size(BZ,2);
obj.f20.Val = [BF(:) BZ(:)];
%EOF
end


