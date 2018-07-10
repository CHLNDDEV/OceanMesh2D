function [] = writefort13( fort13dat, finame )
%
%

if ( nargin == 1 ) 
  finputname = 'fort.13_1' ;
else
  finputname = finame ;   
end

fid = fopen(finputname,'w') ;

% Read a title
% agrid = fgetl(fid) ;
% disp(agrid) ;

fprintf( fid, '%s\n', fort13dat.AGRID ) ; 

fprintf( fid, '%d\n', fort13dat.NumOfNodes ) ; 
fprintf( fid, '%d\n', fort13dat.nAttr ) ;

% Get Default value
nAttr = fort13dat.nAttr ;
for i = 1: nAttr
    % Name
    fprintf( fid, '%s\n', fort13dat.defval.Atr(i).AttrName ) ;
    
    % unit  
    fprintf( fid, '%s\n', fort13dat.defval.Atr(i).Unit ) ;
    
    % values/node 
    fprintf( fid, '%d\n', fort13dat.defval.Atr(i).ValuesPerNode ) ;
    
    % val = fscanf( fid, '%g ', valpernode ) ; 
    
    fprintf( fid, '%15.9e ', fort13dat.defval.Atr(i).Val ) ;
    fprintf( fid, '\n' ) ; 
end

% Get user-defined value
for i = 1: nAttr
    % Name
    fprintf( fid, '%s\n', fort13dat.userval.Atr(i).AttrName ) ;
    
    % No. of Nodes
    fprintf( fid, '%d\n', fort13dat.userval.Atr(i).usernumnodes ) ;
    
    usernodes = fort13dat.userval.Atr(i).usernumnodes ;
    if ( usernodes > 0 )
        % Write user-defined values
        [valpernode,~] = size( fort13dat.userval.Atr(i).Val ) ;
        
        val = fort13dat.userval.Atr(i).Val ;
        
        % Get a proper format
        str = '%d' ;
        for ll = 1: valpernode - 1
            str = [str ' %15.9e'] ;
        end
        str = [str '\n' ] ;
        
        fprintf( fid, str, val ) ;
    end
end

fclose(fid) ; 
