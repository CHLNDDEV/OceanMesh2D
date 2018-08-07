function obj = tidal_data_to_ob(obj,tidal_database,const)
% obj = tidal_data_to_ob(obj,tidal_database,const)
% Input a msh class obj with open boundary locations and tidal
% constituents, and interpolate the tidal database constituents onto it. 
% Put the result into the f15 struct of the msh obj.    
% 
% Requires: m_map for projection                  
%                                                                       
% Created by William Pringle March 15 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check entry
if ~exist(tidal_database,'file')
   error(['tidal database file does not exist: ' tidal_database])        
end
if isempty(obj.op)
   error('No open boundaries in the mesh to put BCs on')
end

% Select desired projection (using m_map)
proj = 'Mercator';
              
% Specify limits of grid
lat_min = min(obj.p(:,2)) - 1; lat_max = max(obj.p(:,2)) + 1;
lon_min = min(obj.p(:,1)) - 1; lon_max = max(obj.p(:,1)) + 1;

% doing the projection
m_proj(proj,'lon',[ lon_min lon_max],...
            'lat',[ lat_min lat_max])

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

%% Load tide data and make vectors
lon = ncread(tidal_database,'lon_z');
lat = ncread(tidal_database,'lat_z');
const_t = ncread(tidal_database,'con');
lon(lon > 180) = lon(lon > 180) - 360;
lon_x = reshape(lon,[],1);
lat_y = reshape(lat,[],1);

% Delete uncessecary portions
% First delete by square
I = find(lon_x < lon_min | lon_x > lon_max | ...
         lat_y < lat_min | lat_y > lat_max);
lon_x(I) = []; lat_y(I) = []; 

% Get 20 nearest neighbours
Kd = ourKNNsearch([lon_x,lat_y]',[b_lon, b_lat]',20);
Kd = unique(Kd);

% The new lon and lat vectors of data
lon_x = lon_x(Kd); lat_y = lat_y(Kd);

% Do the projection
[x,y] = m_ll2xy(lon_x,lat_y);        

%% Now interpolate into f15 struct
keep = true(obj.f15.nbfr,1);
for j = 1:obj.f15.nbfr
    % Read the current consituent
    % For real part
    k = find(startsWith(string(const_t'),lower(const{j})));
    if isempty(k)
       disp(['No tidal data in file for constituent ' const{j}])
       keep(j) = false;  
       continue
    end
    Re_now = ncread(tidal_database,'hRe',[1 1 k],[size(lon) 1]);
    % reshape to vector
    Re_now = reshape(Re_now,[],1);
    % For imaginary part
    Im_now = ncread(tidal_database,'hIm',[1 1 k],[size(lon) 1]);
    % reshape to vector    
    Im_now = reshape(Im_now,[],1);
    % Eliminate regions outside of search area and on land
    % Linear extrapolation of ocean values will be conducted where 
    % boundary nodes fall inside a land cell of the tidal data. 
    Re_now(I) = []; Re_now = Re_now(Kd); 
    K = find(Re_now == 0); Re_now(K) = []; 
    Im_now(I) = []; Im_now = Im_now(Kd); Im_now(K) = [];   
    xx = x; yy = y; xx(K) = []; yy(K) = []; 
    % Make into complex number
    Z = Re_now - Im_now*1i;
    % Do the scattered Interpolation
    F = scatteredInterpolant(xx,yy,Z,'natural');
    BZ = F(b_x,b_y);  
        %
    % Convert real and imaginary parts to amplitude and phase
    amp_b = abs(BZ);  
    phs_b = rad2deg(angle(BZ));
    % Convert to 0 to 360;
    % phs_b that is positive 0 - 180 stays same.
    % phs_b that is negative -180 - 0 becomes 180 - 360
    phs_b(phs_b < 0) = phs_b(phs_b < 0) + 360;
    
    obj.f15.opealpha(j).name = const{j};
    obj.f15.opealpha(j).val = [amp_b phs_b]; 
    
end
obj.f15.nbfr = length(find(keep));
obj.f15.opealpha = obj.f15.opealpha(keep);
obj.f15.bountag = obj.f15.bountag(keep);
%EOF
end
