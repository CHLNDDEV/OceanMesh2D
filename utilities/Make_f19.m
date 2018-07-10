function obj = Make_f19(obj,File_Prefix,File_Suffixes,ETIMINC)
% obj = Make_f19(obj,File_Prefix,File_Suffixes,ETIMINC)
% Input a msh class object get the values of the elevations for times
% between ts and te based on the file_prefix of the input files
%
%  Author:      William Pringle                                 
%  Created:     March 30 2018                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(obj.op)
    error('No boundary information to interpolate to')
end

%% Do some projection for obj.p
proj = 'Mercator';
m_proj(proj,'lon',[ min(obj.p(:,1)) max(obj.p(:,1)) ],...
            'lat',[ min(obj.p(:,2)) max(obj.p(:,2))]) 

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
for t = 1:length(File_Suffixes)
    filename = strcat(File_Prefix,File_Suffixes(t),".nc");
    %% Read the grid and density data 
    if t == 1
        lon = ncread(filename,'lon');
        lat = ncread(filename,'lat');
        lon(lon > 180) = lon(lon > 180) - 360;
        [lon,lat] = ndgrid(lon,lat);
        lon_x = reshape(lon,[],1);
        lat_y = reshape(lat,[],1);
        % Delete uncessecary portions
        % First delete by square
        I = find(lon_x < min(obj.p(:,1)) | lon_x > max(obj.p(:,1)) | ...
                 lat_y < min(obj.p(:,2)) | lat_y > max(obj.p(:,2)));
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
    zeta = ncread(filename,'surf_el');
    Z = reshape(zeta,[],1);
    Z(I) = []; Z = Z(Kd);
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),Z(~isnan(Z)),'natural');
    BZ(:,t) = F(b_x,b_y);  
    
%     scatter(b_lon,b_lat,[],BZ(:,t),'filled')
%     caxis([-0.6 0.3])
%     pause(0.1);
end

%% Make into f19 struct
filename1 = strcat(File_Prefix,File_Suffixes(1),".nc");
filename2 = strcat(File_Prefix,File_Suffixes(end),".nc");
[~,n1] = fileparts(char(filename1)); [~,n2] = fileparts(char(filename2));
obj.f19.DataTitle = [n1 ' -> ' n2];
obj.f19.ETIMINC = ETIMINC;
obj.f19.NumOfNodes = size(BZ,1);
obj.f19.NumOfTimes = size(BZ,2);
obj.f19.Val = BZ(:);
%EOF
end


