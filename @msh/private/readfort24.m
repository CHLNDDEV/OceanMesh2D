function [f24dat] = readfort24( f24name, ncon, nn, salconext )
%
%  f24  -- fort.24 file name
%  salcontext -- sal consituents to be extract (struct)
%
if ( isempty(f24name) )
    f24 = 'fort.24' ;
else
    f24 = f24name;
end

fid = fopen(strtrim(f24)) ;

next = length(salconext) ;

nfound = 0 ;

f24dat.con = salconext ;
f24dat.omega = zeros(next,1) ;
f24dat.found = zeros(next,1) ;
f24dat.samp = zeros(nn,next) ;
f24dat.sphs = zeros(nn,next) ;

for i = 1: ncon
    line = fgetl(fid) ;
    
    cname = strtrim(sscanf(line,'%s %*s')) 
    
    found = 0 ;
    for ii = 1: next
        if ( strcmpi( cname, strtrim(salconext{ii}) )  )
            str = ['Extract ' cname] ; 
            %
            disp(str) ; 
            
            nfound = nfound + 1 ;
            found = 1 ;
            icon = ii ;
            break ;
        end
    end
    
    w = fscanf(fid,'%f \n',1) ;
    line = fgetl(fid) ;
    line = fgetl(fid) ;
    
    val = fscanf(fid,'%d %f %f \n', [3 nn])'  ;
    
    val(1:20,:) 
    if ( found )
        f24dat.omega(icon) = w ;
        f24dat.found(icon) = 1 ;
        f24dat.tiponame{icon} = char(salconext{icon}) ;
        
        f24dat.samp(:,icon) = val(:,2) ;
        f24dat.sphs(:,icon) = val(:,3) ;
    end
    
    if ( nfound == next )
        disp('Done Extracting') ; 
        break ;
    end
    %
end

fclose(fid) ;
return ;

end