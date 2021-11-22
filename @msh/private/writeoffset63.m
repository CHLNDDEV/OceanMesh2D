function fid = writeoffset63( f63dat, finame )
%
%
if ( nargin == 1 ) 
	finputname = 'offset.63' ;
else
	finputname = finame ;   
end

fid = fopen(finputname,'w') ;
fprintf(fid,'%s \n',f63dat.header) ; 
fprintf(fid,'%12.6f \n',f63dat.time_interval) ;
fprintf(fid,'%12.6f \n',f63dat.default_val) ;

for tt = 1: f63dat.num_times
    fprintf(fid,'%d \t %12.6f \n',[f63dat.offset_nodes; f63dat.offset_values(tt,:)]);
    fprintf(fid,'%s \n','##'); 
end
  
fclose(fid) ; 
%EOF
end
