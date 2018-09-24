function outname = writefort14( outfiname, EToV, VX, B, opedat, boudat, title)

% Read ADCIRC grid file fort.14

disp('write adcirc fort.14') ;
tic
%if ( nargin == 0 )
%  finputname = 'fort.14'
%else
%  finputname = finame ;
%end

if ( strcmpi(strtrim(outfiname),'fort.14') )
    disp('output file must not be "fort.14"') ;
    return ;
end

fid = fopen(outfiname,'w') ;
outname = outfiname ;

% agrid = fgetl(fid) ;
% disp(agrid) ;
% title = agrid ;
disp( title )  ;
fprintf(fid,'%s \n', title ) ;

% print NE, NP
NE = length(EToV) ;
NP = length(VX) ;

% write NE, NP
fprintf(fid, '%d %d \n', NE, NP ) ;

% Improve write efficiency
% Vertices
pval = [ [1:1:NP]' VX B] ;
fprintf( fid, '%10d %16.10f %16.10f %18.10e \n', pval' ) ;


% Element connectivity
ElemCon = [ [1:1:NE]' linspace(3,3,NE)' EToV] ;
fprintf(fid, '%d %d %d %d %d \n', ElemCon' ) ;

% Open boundary
if ~isempty(opedat)
    fprintf(fid, '%d %s \n', opedat.nope, '= Number of open boundaries' ) ;
    fprintf(fid, '%d %s \n', opedat.neta, '= Total number of open boundary nodes' ) ;
    
    for i = 1: opedat.nope
        fprintf(fid, ' %d %s \n', opedat.nvdll(i),...
            ['= Number of nodes for open boundary ' num2str(i)]) ;
        
        fprintf( fid, '%d \n', nonzeros(opedat.nbdv(:,i)) ) ;
    end
else
    fprintf(fid, '%d %s \n', 0, '= Number of open boundaries' ) ;
    fprintf(fid, '%d %s \n', 0, '= Total number of open boundary nodes' ) ;
end


% Land boundary
if ~isempty(boudat)
    fprintf(fid, '%d %s \n', boudat.nbou, '= Number of land boundaries' ) ;
    fprintf(fid, '%d %s \n', boudat.nvel, '= Total number of land boundary nodes' ) ;
    
    %
    for i = 1: boudat.nbou
        % Bc type:q
        fprintf( fid, '%d %d %s \n', boudat.nvell(i), boudat.ibtype(i), ...
            ['= Number of nodes for land boundary ' num2str(i)]) ;
        
        switch ( boudat.ibtype(i) )
            case {0,1,2,10,11,12,20,21,22,30,60,61,101,52}
                % for k = 1: nvell(i)
                %   msgline  = fgetl(fid) ;
                %   nbvv(k,i) = str2num(msgline) ;
                % end
                % Nov 15, 2012, improve reading efficiency
                % nbvv(1:nvell(i),i) = fscanf(fid,'%d \n', nvell(i) ) ;
                fprintf( fid, '%d \n' , nonzeros(boudat.nbvv(:,i)) ) ;
            case  {3, 13, 23}
                disp('Warning: 3 13 23 --- To be added')
                % val = fscanf(fid,'%g %g %g \n', [3 nvell(i)] )  ;
                % nbvv(1:nvell(i),i) = val(1,:) ;
            case  {4, 24}
                % disp('4 24')
                % val = fscanf(fid,'%g %g %g %g %g \n', [5 nvell(i)] )  ;
                idx = find( boudat.nbvv(:,i) ) ;
                
                nbvv = full(boudat.nbvv(idx,i)) ;
                ibconn = full(boudat.ibconn(idx,i)) ;
                barinht = full(boudat.barinht(idx,i)) ;
                barincfsb = full(boudat.barincfsb(idx,i)) ;
                barincfsp = full(boudat.barincfsp(idx,i)) ;
                for ll = 1: length(idx)
                    %
                    % fprintf( fid, '%d %d %f %f %f \n', boudat.nbvv(idxn,i), ...
                    %    boudat.ibconn(idxn,i), boudat.barinht(idxn,i), ...
                    %    boudat.barincfsb(idxn,i), boudat.barincfsp(idxn,i) ) ;
                    fprintf( fid, '%d %d %16.10f %16.10f %16.10f \n', nbvv(ll), ...
                        ibconn(ll), barinht(ll), ...
                        barincfsb(ll), barincfsp(ll) ) ;
                    %
                end
                
            case  {5, 25}
                disp('Warning: 5 25 --- to be added')
                % val = fscanf(fid,'%g % g %g %g %g %g %g %g \n', [8 nvell(i)] ) ;
                % nbvv(1:nvell(i),i) = val(1,:) ;
                %otherwise
                %    msgline = fgetl(fid) ;
            case 94
                fprintf( fid, '%d %d \n' , boudat.nbvv(:,1:2)' ) ;
                
        end
    end
    % case of no bou
else
    fprintf(fid,'%d %s\n',0,' !');
    fprintf(fid,'%d %s\n',0,' !');
    fprintf(fid,'%d %s\n',0,' !');
    fprintf(fid,'%d %s\n',0,' !');
    
end
% toc

fclose(fid) ;

return

