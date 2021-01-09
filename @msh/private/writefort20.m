function [] = writefort20( fort20dat, finame )
%
%

if ( nargin == 1 ) 
  finputname = 'fort.20_1' ;
else
  finputname = finame ;   
end

fid = fopen(finputname,'w') ;

fprintf( fid, '%d\n', fort20dat.FTIMINC ) ; 

valpernode = size( fort20dat.Val,2 ) ;
str = [];
for ll = 1: valpernode
    str = [str '%.8f '] ;
end
str = [str '\n' ] ;

fprintf( fid, str, fort20dat.Val' ) ;
    
fclose(fid) ; 
%EOF
end
