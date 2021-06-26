function fid = writefort24( f24dat, finame, format )
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
    % SAL indices
    indices = squeeze(f24dat.Val(1,1,:)); 
    sal_nodes = length(indices);
    nccreate(file,'mesh_indices','Dimensions',{'sal_nodes',sal_nodes},'DataType','int32','DeflateLevel',dl) ;
    ncwriteatt(file,'mesh_indices','long_name','index of SAL values into mesh coordinates')
    ncwrite(file, 'mesh_indices', indices) ;
    % SAL name
    tipnames = char(f24dat.tiponame); ntip = size(tipnames) ;
    nccreate(file,'const','dimensions',{'num_const',ntip(1),'char_len',ntip(2)},'DataType','char') ;
    ncwriteatt(file,'const','long_name','names of the tidal harmonic constituents')
    ncwrite(file, 'const', tipnames) ;
    % SAL frequencies
    nccreate(file,'frequency','dimensions',{'num_const',ntip(1)});
    ncwriteatt(file,'frequency','long_name','frequency of harmonic constituents')
    ncwriteatt(file,'frequency','units','rad/s')
    ncwrite(file, 'frequency', f24dat.omega) ;
    % SAL amplitudes
    nccreate(file, 'salamp','Dimensions',{'num_const',ntip(1),'sal_nodes',sal_nodes},'DataType','single','DeflateLevel',dl) ;
    ncwriteatt(file,'salamp','long_name','amplitude of self-attraction and loading tide elevation')
    ncwriteatt(file,'salamp','units','m')
    ncwrite(file, 'salamp', squeeze(f24dat.Val(:,2,:))) ;
    % SAL phase lags
    nccreate(file, 'salphs','Dimensions',{'num_const',ntip(1),'sal_nodes',sal_nodes},'DataType','single','DeflateLevel',dl) ;
    ncwriteatt(file,'salphs','long_name','phase-lag of self-attraction and loading tide elevation')
    ncwriteatt(file,'salphs','units','degrees with respect to GMT/UTC')
    ncwrite(file, 'salphs', squeeze(f24dat.Val(:,3,:))) ;
    % Add some meta data/attribution
    ncwriteatt(file,'/','title','The self-attraction and loading terms for an ADCIRC simulation');
    ncwriteatt(file,'/','creation_date',datestr(now));
    ncwriteatt(file,'/','source',"Made by OceanMesh2D writefort24");
    ncwriteatt(file,'/','references',"https://github.com/CHLNDDEV/OceanMesh2D/" );

else
    error(['format = ' format ' is invalid. Choose from ascii or netcdf'])
end
%EOF
end
