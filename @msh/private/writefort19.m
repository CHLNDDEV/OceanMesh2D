function [] = writefort19( fort19dat, finame )
%
%

if ( nargin == 1 ) 
  finputname = 'fort.19_1' ;
else
  finputname = finame ;   
end

fid = fopen(finputname,'w') ;

fprintf( fid, '%d\n', fort19dat.ETIMINC ) ; 

valpernode = size( fort19dat.Val,2 ) ;
str = [];
for ll = 1: valpernode
    str = [str '%.8f '] ;
end
str = [str '\n' ] ;

fprintf( fid, str, fort19dat.Val' ) ;
    
fclose(fid) ; 
%EOF
end
