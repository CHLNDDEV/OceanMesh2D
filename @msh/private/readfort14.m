function [EToV,VX,B,opedat,boudat,title] = readfort14( finame, read_bou )
% Read ADCIRC grid file fort.14

disp('read adcirc fort.14') ;
tic
if ( nargin == 0 )
    finputname = 'fort.14'
else
    finputname = finame ;
end

% open file
fid = fopen(finputname) ;

% read first title line
agrid = fgetl(fid) ;
disp(agrid) ;
title = agrid ;

% read number of nodes and elements
msgline = fgetl(fid) ;
N = sscanf(msgline,'%d %d %*s') ;

% read the node numbers, position and depth
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;

% read the triangulation
idx = fscanf(fid,'%d %d %d %d %d \n', [5 N(1)])' ;

% sort read items into vectors/matrices
EToV = idx(:,3:5) ;
% treat non-sequential ordering
VX = NaN(max(EToV(:)),2);
B = NaN(max(EToV(:)),1);
VX(Val(:,1),:) = Val(:,2:3) ;
B(Val(:,1)) = Val(:,4) ;

if read_bou
    % Read in boundary
    % Open boundary
    msgline = fgetl(fid) ;
    nope = sscanf(msgline,'%d %*s') ;
    
    msgline = fgetl(fid) ;
    neta = sscanf(msgline,'%d %*s') ;
    
    nvdll = zeros(nope,1) ;
    ibtypee = zeros(nope,1) ;
    nbdv = zeros(neta,nope) ;
    % tic
    for i = 1: nope
        msgline = fgetl(fid) ;
        
        [varg] = sscanf(msgline,'%d %*s \n') ;
        nvdll(i) = varg ;
        ibtypee(i) = 0 ;
        
        %
        % for k = 1: nvdll(i)
        %
        %    % % nbdv(k,i) = fscanf(fid,'%d \n')
        %    % % fscanf(fid,'%d \n')
        %    msgline = fgetl(fid) ;
        %
        %    nbdv(k,i) = str2num(msgline) ;
        % end
        %
        % Nov 25, 2012, improve reading efficiency
        nbdv(1:nvdll(i),i) = fscanf(fid,'%d \n', nvdll(i) ) ;
    end
    % toc
    
    % ocean boundary
    opedat.nope = nope ;
    opedat.neta = neta ;
    opedat.nvdll = nvdll ;
    opedat.ibtypee = ibtypee ;
    opedat.nbdv = nbdv ;
    
    % land boundary
    msgline = fgetl(fid) ;
    nbou = sscanf(msgline,'%d %*s') ;
    
    msgline = fgetl(fid) ;
    nvel = sscanf(msgline,'%d %*s') ;
    
    nvell = zeros(nbou,1) ;
    ibtype = zeros(nbou,1) ;
    nbvv = sparse(nvel,nbou) ;
    ibconn = sparse(nvel,nbou) ;
    barinht = sparse(nvel,nbou) ;
    barincfsb = sparse(nvel,nbou) ;
    barincfsp = sparse(nvel,nbou) ;
    
    % tic
    for i = 1: nbou
        msgline = fgetl(fid) ;
        
        [varg,~] = sscanf(msgline,'%d %d %*s \n') ;
        nvell(i) = varg(1) ;
        ibtype(i) = varg(2) ;
        
        switch ( ibtype(i) )
            case {0,1,2,10,11,12,20,21,22,30,60,61,101,52}
                % Nov 15, 2012, improve reading efficiency
                nbvv(1:nvell(i),i) = fscanf(fid,'%d \n', nvell(i) ) ;
            case  {3, 13, 23}
                disp('3 13 23')
                val = fscanf(fid,'%g %g %g \n', [3 nvell(i)] )  ;
                nbvv(1:nvell(i),i) = val(1,:) ;
            case  {4, 24}
                %disp('4 24')
                val = fscanf(fid,'%g %g %g %g %g \n', [5 nvell(i)] )  ;
                nbvv(1:nvell(i),i) = val(1,:) ;
                ibconn(1:nvell(i),i) = val(2,:) ;
                barinht(1:nvell(i),i) = val(3,:) ;
                barincfsb(1:nvell(i),i) = val(4,:) ;
                barincfsp(1:nvell(i),i) = val(5,:) ;
            case  {5, 25}
                %disp('5 25')
                val = fscanf(fid,'%g % g %g %g %g %g %g %g \n', [8 nvell(i)] ) ;
                nbvv(1:nvell(i),i) = val(1,:) ;
                %otherwise
                %    msgline = fgetl(fid) ;
            case  94
                val = fscanf(fid,'%d %d \n', [2 nvell(i)] ) ;
                nbvv(1:nvell(i),1:2) = val' ;
        end
    end
    % toc
    
    % land boundary
    boudat.nbou = nbou ;
    boudat.nvel = nvel ;
    boudat.nvell = nvell ;
    boudat.ibtype = ibtype ;
    boudat.nbvv = nbvv ;
    
    if ( sum(ibtype == 24) > 0 ||  sum(ibtype == 4) > 0 )
        boudat.ibconn = ibconn ;
        boudat.barinht = barinht ;
        boudat.barincfsb = barincfsb ;
        boudat.barincfsp = barincfsp ;
    end
else
    opedat = []; 
    boudat  = []; 
end

fclose(fid) ;

return

