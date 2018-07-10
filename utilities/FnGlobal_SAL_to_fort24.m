function fid = FnGlobal_SAL_to_fort24( f24out, grd, avisoloc, saldata )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolates the global SAL term onto the mesh and outputs a fort.24    %                                                                       %
%                                                                         %  
% Requires: readfort14.m                                                  %
%                                                                         % 
% Data required:                                                          %
% FES2004 loads. Source at: ftp://ftp.legos.obs-mip.fr/pub/soa/maree/...  %
%                           tide_model/global_solution/fes2004/           %
%                                                                         %
% FES2014 loads. Source at: ftp://ftp.legos.obs-mip.fr/pub/FES2012-project/data/LSA/FES2014/             %
%                                                                         %
% Created by William Pringle Oct 20 2016 for FES2004 SAL                  %
% Updated by William Pringle Oct 28 2016 for FES2014 SAL                  %
% Updated by Dam Wa for make into a function                              %
%                                                                         %
% Run example:                                                            %
% f14 = 'fort.14' ;                                                       %
% f24 = 'fort.24'                                                         %
% latlon0 = [75.214667 -31.172085] ;                                      %
% saldat = 'FES2014' ;                                                    %
% avisoloc = './AVISO_DIREC'                                              %
% f15tipname = { 'M2', 'O1', 'S2', 'N2', 'K2', 'K1', 'Q1', 'P1'}          %
%                                                                         %
% Global_SAL_to_fort24( f24, f14, f15tipname, lonlat0, avisoloc, saldata )%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll0 = grd.f15.slam(1) ;
if ( ll0 < 0 ) 
    ll0 = ll0 + 360 ; 
end
lon0 = ll0*pi/180 ; lat0 = grd.f15.slam(2)*pi/180 ; 

R = 6378206.4; % earth radius  

f15tipname = {grd.f15.tipotag.name};

ntip = length(f15tipname) ;  

% choose tidal database file names and directories
%database = 'FES2004';
% database = 'FES2014';
database = strtrim(upper(saldata)) ;
% direc    = 'E:\Global_Data\AVISO_TIDES\';
direc = strtrim(avisoloc) ;

% output fort.24 name
% fort24    = ['fort.24.' database];
fort24 = [strtrim(f24out) '.' database] ;

% % Load tide grid data 
if strcmp(database,'FES2004')
    tide_grid     = [direc '/' database '/SAL/load.k1.nc'];
    tide_prefix   = [direc '/' database '/SAL/load.'];
    tide_suffix   = '.nc';
    if ( ispc )
        tide_grid     = [direc '\' database '\SAL\load.k1.nc'];
        tide_prefix   = [direc '\' database '\SAL\load.'];
        % tide_suffix   = '.nc';
    end
    
    %ncdisp(tide_grid);
    lon = ncread(tide_grid,'lon');
    lat = ncread(tide_grid,'lat');
elseif  strcmp(database,'FES2014')
    tide_grid     = [direc '/' database '/SAL/K1_sal.nc'];
    tide_prefix   = [direc '/' database '/SAL/'];
    tide_suffix   = '_sal.nc';
    if ( ispc )
        tide_grid     = [direc '\' database '\SAL\K1_sal.nc'];
        tide_prefix   = [direc '\' database '\SAL\'];
        tide_suffix   = '_sal.nc';
    end
    
    %ncdisp(tide_grid);
    lon = ncread(tide_grid,'longitude');
    lat = ncread(tide_grid,'latitude');
    [lon,lat] = ndgrid(lon,flipud(lat));
end
% CPP Conversion of lat/lon
lon = lon * pi/180; lat = lat * pi/180; 
x = R * (lon - lon0) * cos(lat0);
y = R * lat;

% % Get mesh info
VX = grd.p;

% Doing the CPP conversion
VX(VX(:,1)<0,1)=VX(VX(:,1)<0,1)+360; 

xx = VX(:,1) * pi/180; yy = VX(:,2) * pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

% % Now interpolate onto grid and write out to fort.24 type file

fid = fopen(fort24,'w');

nnodes = length(VX) ;
kvec = [1:nnodes]' ; 
for icon = 1: ntip
% for j = 1:length(const)

    % The current consituent filename
    if strcmp(database,'FES2004')
        tide = [tide_prefix lower(f15tipname{icon}) tide_suffix];
        % Get amp and phase
        Ha = ncread(tide,'Ha');
        Hg = ncread(tide,'Hg');
    elseif  strcmp(database,'FES2014')
        tide = [tide_prefix f15tipname{icon} tide_suffix];
        % Get amp and phase
        Ha = ncread(tide,'SAL_amplitude');
        Hg = ncread(tide,'SAL_phase');
        Ha = fliplr(Ha);
        Hg = fliplr(Hg);
    end

%     figure;
%     [c,h] = contour(lon,lat,Ha);
%     clabel(c,h)
%     colorbar
%     figure;
%     [c,h] = contour(lon,lat,Hg,0:30:360);
%     clabel(c,h)
%     colorbar
    
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
    
    % Print out interpolated results
    
    fprintf('Writing SAL %s data \n', char(f15tipname{icon})) ; 
    % The constituent details
    fprintf(fid,'%s \n',[char(f15tipname{icon}) ' SAL']) ;
    fprintf(fid,'%17.15f \n',grd.f15.tipotag(icon).val(2)) ;
    fprintf(fid,'%d \n',1) ;  
    fprintf(fid,'%s \n',char(f15tipname{icon})) ;
    
    % Loop over the nodes of the mesh
    % for k = 1:length(amp)
    %    fprintf(fid,'%d \t %12.6f  %12.6f \n',k,amp(k),phs(k));
    % end
    figure; fastscatter(xx,yy,amp); 
    pause(0.5)
    fprintf(fid,'%d \t %12.6f  %12.6f \n',[kvec amp phs]');
    
end
fclose(fid);
end
