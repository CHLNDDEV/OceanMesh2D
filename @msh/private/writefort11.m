function [] = writefort11( fort11dat, finame )
%
%

if ( nargin == 1 ) 
  finputname = 'fort.11_1' ;
else
  finputname = finame ;   
end

fid = fopen(finputname,'w') ;

fprintf( fid, '%s\n', fort11dat.DataTitle ) ; 
%fprintf( fid, '%s\n', fort11dat.DataSubTitle ) ; 
fprintf( fid, '%d\n', fort11dat.DTIMINC ) ; 

fprintf( fid, '%d\n', fort11dat.NumOfNodes ) ; 

if iscell(fort11dat.Val)
    ValL = length(fort11dat.Val);
    valpernode = size(fort11dat.Val{1},1);
else
    % Write user-defined values
    [valpernode,~,ValL] = size( fort11dat.Val ) ;
end

% Looping over the number of times
for t = 1:ValL

    if iscell(fort11dat.Val)
        val = fort11dat.Val{t};
    else
        val = fort11dat.Val(:,:,t) ;
    end
    % Get a proper format
    str = '%d' ;
    for ll = 1: valpernode - 1
        str = [str ' %15.9e'] ;
    end
    str = [str '\n' ] ;

    fprintf( fid, str, val ) ;
end
    
fclose(fid) ; 
%EOF
end
