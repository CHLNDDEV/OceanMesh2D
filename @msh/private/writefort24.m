function fid = writefort24( f24dat, finame, ascii )
%
%
if ( nargin == 1 )
    finputname = 'fort.24_1' ;
else
    finputname = finame ;
end

if nargin < 3
    ascii = 0;
end


if ascii % LEGACY
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
    
else
    % create nc file or overwrite it
    file = [finputname,'.nc'];
    if exist(file)
        delete(file)
    end
    nnodes = length(f24dat.Val(1,:,:));
    tipnames = f24dat.tiponame; ntip = length(tipnames) ;
    for icon = 1: ntip
        tipname = tipnames{icon};
        var1= [tipname,'_omega'];
        fprintf('Writing SAL %s data \n', char(tipname)) ;
        nccreate(file,var1) ;
        ncwrite(file, var1, f24dat.omega(icon)) ;
        
        var2= [tipname,'_vals'];
        nccreate(file, var2,'Dimensions',{'valspernode',3,'nnodes',nnodes}) ;
        ncwrite(file, var2, squeeze(f24dat.Val(icon,:,:)),[1, 1]) ;
    end
    
end
%EOF
end