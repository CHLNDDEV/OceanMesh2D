function fid = writefort24( f24dat, finame )
%
%
if ( nargin == 1 ) 
	finputname = 'fort.24_1' ;
else
	finputname = finame ;   
end

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
%EOF
end