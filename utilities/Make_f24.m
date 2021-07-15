function obj = Make_f24( obj, saldata, plot_on )
% obj = Make_f24( obj, saldata, plot_on )
% Takes a msh object and interpolates the global SAL term into the f24
% struct
% Assumes that saldata is in the MATLAB path
% The saldata required can be downloaded from:
% saldata = 'FES2004' : Source at: ftp://ftp.legos.obs-mip.fr/pub/soa/...
%                                 maree/tide_model/global_solution/fes2004/load/
%
% saldata = 'FES2014' : Source at: ftp://ftp.legos.obs-mip.fr/pub/...
%                                  FES2012-project/data/LSA/FES2014/
% by default saldata = 'FES2014'
%
% plot_on -  1/true: to plot and print F24 values for checking
%            0/false: no plotting by default
%
% Created by William Pringle. July 11 2018 updated to Make_f## style
% Modified by Keith Roberts Jun 19, 2021 to specify degree format for FES
% file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(obj.f15)
    error(['The msh object must have the f15 struct populated '...
        'with tidal potential information'])
end

if nargin < 2 || isempty(saldata)
    saldata = 'FES2014';
end
if nargin < 3 || isempty(plot_on)
    plot_on = false;
end

% parameter for cpp conversion
R = 6378206.4; % earth radius
lon0 = obj.f15.slam(1) ; % central longitude
lat0 = obj.f15.slam(2) ; % central latitude
% Put in the tidal potential names
obj.f24.tiponame = {obj.f15.tipotag.name};

ntip = length(obj.f24.tiponame) ;

% choose tidal database file names and directories
database = strtrim(upper(saldata)) ;
direc    = '';

% % Load tide grid data
if strcmp(database,'FES2004')
    tide_grid     = [direc 'load.k1.nc'];
    tide_prefix   = [direc 'load.'];
    tide_suffix   = '.nc';
    lon = ncread(tide_grid,'lon');
    lat = ncread(tide_grid,'lat');
    % -180/180 degree format for 2004
    if (lon0 > 180); lon0 = lon0 - 360 ; end
elseif strcmp(database,'FES2014')
    tide_grid     = [direc  'K1_sal.nc'];
    tide_prefix   = direc;
    tide_suffix   = '_sal.nc';
    lon = ncread(tide_grid,'longitude');
    lat = ncread(tide_grid,'latitude');
    [lon,lat] = ndgrid(lon,flipud(lat));
    % 0-360 degree format for 2014
    if (lon0 < 0); lon0 = lon0 + 360 ; end
else
    error(['database = ' database ' is invalid.'...
          ' Choose from FES2004 or FES2014'])
end
% Convert CPP origin parameters to radians
lon0 = lon0*pi/180 ; lat0 = lat0*pi/180 ;

% CPP Conversion of FES lat/lon
lon = lon * pi/180; lat = lat * pi/180;
x = R * (lon - lon0) * cos(lat0);
y = R * lat;

% Doing the CPP conversion of the mesh
VX = obj.p;
if strcmp(database,'FES2004')
   VX(VX(:,1)>180,1) = VX(VX(:,1)>180,1) - 360;
elseif strcmp(database,'FES2014')
   VX(VX(:,1)<0,1) = VX(VX(:,1)<0,1) + 360;
end
xx = VX(:,1) * pi/180; yy = VX(:,2) * pi/180;
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

% Now interpolate onto grid and put into fort.24 struct
nnodes = length(VX) ;
kvec = (1:nnodes)';
obj.f24.Val = zeros(ntip,3,nnodes) ;
for icon = 1: ntip
    % Get the frequency
    obj.f24.omega(icon) = obj.f15.tipotag(icon).val(2);
    % The current consituent filename
    if strcmp(database,'FES2004')
        tide = [tide_prefix lower(obj.f24.tiponame{icon}) tide_suffix];
        % Get amp and phase
        Ha = ncread(tide,'Ha');
        Hg = ncread(tide,'Hg');
    elseif  strcmp(database,'FES2014')
        tide = [tide_prefix obj.f24.tiponame{icon} tide_suffix];
        % Get amp and phase
        Ha = ncread(tide,'SAL_amplitude');
        Hg = ncread(tide,'SAL_phase');
        Ha = fliplr(Ha);
        Hg = fliplr(Hg);
    end
    
    Hg(Hg > 180) = Hg(Hg > 180) - 360; % move to -180 - 180
    Hg = Hg*pi/180; %radians
    
    % Convert to complex number for interpolation
    z = Ha.*exp(Hg*1i);
    
    % Do the gridded Interpolation
    F = griddedInterpolant(x,y,z,'linear','none');
    Z = F(xx,yy);
    
    % Convert back to amp and phase
    amp = abs(Z);
    phs = angle(Z)*180/pi;
    
    % Convert to 0 to 360;
    % phs_b that is positive 0 - 180 stays same.
    % phs_b that is negative -180 - 0 becomes 180 - 360
    phs(phs < 0) = phs(phs < 0) + 360;
    
    % Plot interpolated results
    if plot_on
        figure(1); fastscatter(VX(:,1),VX(:,2),amp);
        colorbar;
        constituent = obj.f24.tiponame{icon};
        title(constituent)
        print(['F24_' constituent '_check'],'-dpng')
    end
    
    % Put into the struct
    obj.f24.Val(icon,:,:) = [kvec'; amp'; phs'];
    
    if any(isnan(amp))
        warning('NaNs detected in amplitudes. Is degree format correct?')
    end
end

if obj.f15.ntip ~= 2
    disp('Setting ntip = 2')
    obj.f15.ntip = 2;
end
%EOF
end
