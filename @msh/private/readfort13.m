function [fort13dat] = readfort13( finame )
%
% 
% Read fort.13
%
%
% DW:
%

if ( nargin == 0 ) 
  finputname = 'fort.13' ;
else
  finputname = finame ;   
end

fid = fopen(finputname) ;

%
% Read a title
agrid = fgetl(fid) ;
disp(agrid) ;

fort13dat.AGRID = agrid ;

%
% Read a number of node
nnodes = fscanf(fid,'%d \n',1) ;
nAttr = fscanf(fid,'%d \n', 1) ; 

fort13dat.NumOfNodes = nnodes ; 
fort13dat.nAttr = nAttr ;

%
% Get Default value
for i = 1: nAttr
    attr = fgetl(fid) ;
    fort13dat.defval.Atr(i).AttrName = strtrim(attr) ;
    
    unit = fgetl(fid) ; 
    fort13dat.defval.Atr(i).Unit = strtrim(unit) ;
    
    valpernode = fscanf(fid, '%d \n', 1 ) ;
    fort13dat.defval.Atr(i).ValuesPerNode = valpernode ;
    
    val = fscanf( fid, '%g ', valpernode ) ;  
    fort13dat.defval.Atr(i).Val = val ;
    v1 = fscanf( fid, '\n' ) ; % flush
end

%
% Get user-defined value
for i = 1: nAttr
% for i = 1: 1
    attr = fgetl(fid) ;
    fort13dat.userval.Atr(i).AttrName = strtrim(attr) ;
    
    numnodes = fscanf(fid, '%d \n', 1 ) ; 
    fort13dat.userval.Atr(i).usernumnodes = numnodes ;
    
    for j = 1: nAttr
    % for j = 1: 1
        fnd = strcmpi(strtrim(fort13dat.defval.Atr(j).AttrName),strtrim(attr)) ;
       
        if ( fnd )
            valpernode = fort13dat.defval.Atr(j).ValuesPerNode ;
            break ; 
        end
    end
    
    %
    if ( fnd ) 
      % Record value
      val = fscanf(fid, '%g ', (valpernode + 1)*numnodes ) ;
      fort13dat.userval.Atr(i).Val = reshape(val,valpernode + 1,numnodes ) ; 

      v1 = fscanf( fid, '\n' ) ; % flush a new line
    else
      errstr = sprintf('Error: attribute %s is not found!', attr ) ;    
      disp(errstr) ; 
    end
end

%
fclose(fid) ;