function fid = writefort24( f24dat, finame, points, format )
%
%
if ( nargin == 1 )
    finputname = 'fort.24_1' ;
else
    finputname = finame ;
end

if nargin < 3 || isempty(format) 
   format = 'ascii';
   warning('Using ASCII file format. Pass ''netcdf'' to write in NetCDF')
end
if iscell(format)
   format = format{1};
end

if strcmp(format,'ascii') % LEGACY
    fid = fopen(finputname,'w') ;
    
    tipnames = f24dat.tiponame; ntip = length(tipnames) ;
    for icon = 1: ntip
        tipname = tipnames{icon};
        fprintf('Writing SAL %s data \n', char(tipname)) ;
        % The constituent details
        fprintf(fid,'%s \n',[char(tipname) ' SAL']) ;
        fprintf(fid,'%17.15f \n',f24dat.omega(icon)) ;
        fprintf(fid,'%d \n',1) ;
        fprintf(fid,'%s \n',char(tipname)) ;
        fprintf(fid,'%d \t %12.6f  %12.6f \n',f24dat.Val(icon,:,:));
    end
    
    fclose(fid) ;
    
elseif strcmp(format,'netcdf') % NETCDF
    % create nc file or overwrite it
    file = [finputname,'.nc'];
    if exist(file)
        delete(file)
    end
    % deflate level (set zero for no deflation if necessary.. using single precision so not very different in size)
    dl = 5;
    
    node = length(points); 
    tipnames = char(f24dat.tiponame); ntip = size(tipnames) ;
    fillvalue = 0.0;
        
    nccreate(file,'x','Dimensions',{'node',node},'DataType','int32','DeflateLevel',dl) ;
    ncwriteatt(file,'x','standard_name', 'longitude') ;
    ncwriteatt(file,'x','units', 'degrees_east') ;
    ncwriteatt(file,'x','positive', 'east') ;
    ncwrite(file, 'x', points(:,1)) ;

    nccreate(file,'y','Dimensions',{'node',node},'DataType','int32','DeflateLevel',dl) ;
    ncwriteatt(file,'y','standard_name', 'latitude') ;
    ncwriteatt(file,'y','units', 'degrees_north') ;
    ncwriteatt(file,'y','positive', 'north') ;
    ncwrite(file, 'y', points(:,2)) ;    
    
    nccreate(file,'constituents','Dimensions',{'num_constituents',ntip(1),'char_len',ntip(2)},'DataType','char') ;
    ncwriteatt(file,'constituents','standard_name', 'name_of_tidal_harmonic_constituents') ;
    ncwriteatt(file,'constituents','long_name', 'name of tidal harmonic constituents') ;
    ncwrite(file, 'constituents', tipnames) ;
    
    nccreate(file,'frequency','dimensions',{'num_constituents',ntip(1)});
    ncwriteatt(file,'frequency','standard_name','frequency_of_tidal_harmonic_constituents')
    ncwriteatt(file,'frequency','long_name','frequency of tidal harmonic constituents')
    ncwriteatt(file,'frequency','units','radians/second')
    ncwrite(file, 'frequency', f24dat.omega) ;
    
    nccreate(file, 'sal_amplitude','Dimensions',{'num_constituents',ntip(1),'node',node},'DataType','single','DeflateLevel',dl,'FillValue',fillvalue) ;
    ncwriteatt(file,'sal_amplitude','standard_name','amplitude_of_self_attraction_and_loading_tide_elevation')
    ncwriteatt(file,'sal_amplitude','long_name','amplitude of self attraction and loading tide elevation')
    ncwriteatt(file,'sal_amplitude','units','m')
    ncwrite(file, 'sal_amplitude', squeeze(f24dat.Val(:,2,:))) ;
    
    nccreate(file, 'sal_phase','Dimensions',{'num_constituents',ntip(1),'node',node},'DataType','single','DeflateLevel',dl,'FillValue',fillvalue) ;
    ncwriteatt(file,'sal_phase','long_name','phase-lag of self-attraction and loading tide elevation')
    ncwriteatt(file,'sal_phase','standard_name','phase_lag_of_self_attraction_and_loading_tide_elevation')
    ncwriteatt(file,'sal_phase','units','degrees with respect to GMT/UTC')
    ncwrite(file, 'sal_phase', squeeze(f24dat.Val(:,3,:))) ;
   
    ncwriteatt(file,'/','title','The self-attraction and loading terms for an ADCIRC simulation');
    ncwriteatt(file,'/','creation_date',datestr(now));
    ncwriteatt(file,'/','source',"Made by OceanMesh2D writefort24");
    ncwriteatt(file,'/','references',"https://github.com/CHLNDDEV/OceanMesh2D/" );

else
    error(['format = ' format ' is invalid. Choose from ascii or netcdf'])
end
%EOF
end
