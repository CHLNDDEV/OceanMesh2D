function [fid] = writefort5354( f5354dat, finame )
%
fort53 = 'fort.53001' ;
fort54 = 'fort.54001' ;
if ( nargin >= 1 )
    fort53 = [strtrim(finame) '.53001'];
    fort54 = [strtrim(finame) '.54001'];
end

% open the files
fid = fopen(fort53, 'w') ;
fid1 = fopen(fort54, 'w') ;

% number of frequency
fprintf(fid, '%d \n', f5354dat.nfreq ) ;
fprintf(fid1, '%d \n', f5354dat.nfreq ) ;

nfreq = f5354dat.nfreq ;
for i = 1: nfreq
    fprintf( fid, '%16.9e %f %f %s \n', f5354dat.freqinfo(i).val, ...
                                    f5354dat.freqinfo(i).name ) ;
    fprintf( fid1, '%16.9e %f %f %s \n', f5354dat.freqinfo(i).val, ...
                                    f5354dat.freqinfo(i).name ) ;
end

% Number of nodes/stations
nstae = length(f5354dat.nodes);
fprintf(fid, '%d \n', nstae) ;
fprintf(fid1,'%d \n', nstae) ;

% write out each node/stations
for i = 1: nstae
    fprintf(fid,'%d \n', f5354dat.nodes(i)) ;  
    fprintf(fid1,'%d \n', f5354dat.nodes(i)) ;  
    
    val = f5354dat.ele(:,:,i) ;
    fprintf(fid,'  %f %f \n', val' ) ;
    
    val = f5354dat.vel(:,:,i) ;
    fprintf(fid1,' %f %f %f %f \n', val' ) ;
end

% close the file
fclose(fid) ; fclose(fid1) ;