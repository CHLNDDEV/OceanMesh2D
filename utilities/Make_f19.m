function obj = Make_f19(obj,filename)
% obj = Make_f19(obj,filename)
% Input a msh class object get the values of the elevations for all times
% in the filename
%
%  Author:      William Pringle                                 
%  Created:     March 30 2018
%               July 19 2018, Updated to do it just for one file                                  
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

time = ncread(filename,'time');
DT = median(diff(time))*3600; %[s]
L = length(time);
BZ = zeros(length(b_x),L);
disp(['Interpolating ' num2str(L) ' snaps'])
for t = 1:L
    if mod(t,100) == 0 
        disp(['Snap #' num2str(t)])
    end
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
        Kd = ourKNNsearch([lon_x,lat_y]',[b_lon, b_lat]',20);
        Kd = unique(Kd);

        % The new lon and lat vectors of data
        lon_x = lon_x(Kd); lat_y = lat_y(Kd);

        % Do the projection
        [x,y] = m_ll2xy(lon_x,lat_y);       
    end
    zeta = ncread(filename,'surf_el',[1 1 t],[size(lon) 1]);
    Z = reshape(zeta,[],1);
    Z(I) = []; Z = Z(Kd);
    F = scatteredInterpolant(x(~isnan(Z)),y(~isnan(Z)),Z(~isnan(Z)),'natural');
    BZ(:,t) = F(b_x,b_y);  
end

%% Make into f19 struct
ts = datenum('2000-01-01') + time(1)/24;
te = datenum('2000-01-01') + time(end)/24;
obj.f19.DataTitle = [datestr(ts) ' -> ' datestr(te)];
obj.f19.ETIMINC = DT;
obj.f19.NumOfNodes = size(BZ,1);
obj.f19.NumOfTimes = size(BZ,2);
obj.f19.Val = BZ(:);
%EOF
end


