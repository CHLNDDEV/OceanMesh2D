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
    indices = squeeze(f24dat.Val(1,1,:)); 
    nnodes = length(indices);
    nccreate(file,'mesh_indices','Dimensions',{'nnodes',nnodes},'DataType','int32','DeflateLevel',7) ;
    ncwrite(file, 'mesh_indices', indices) ;

    tipnames = f24dat.tiponame; ntip = length(tipnames) ;
    for icon = 1: ntip
        tipname = tipnames{icon};
        var1= [tipname,'_omega'];
        fprintf('Writing SAL %s data \n', char(tipname)) ;
        nccreate(file,var1) ;
        ncwrite(file, var1, f24dat.omega(icon)) ;
        
        var2= [tipname,'_vals'];
        nccreate(file, var2,'Dimensions',{'valspernode',2,'nnodes',nnodes},'DeflateLevel',7) ;
        ncwrite(file, var2, squeeze(f24dat.Val(icon,2:3,:)),[1,1]) ;
    end

else
    error(['format = ' format ' is invalid. Choose from ascii or netcdf'])
end
%EOF
end
