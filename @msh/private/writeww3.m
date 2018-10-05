function outname = writeww3(outfiname,EToV,VX,B,opedat, title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes WW3 grid including nodes             %
% (longitude,latitude,depth), open bounday nodes            %
% and element connections (triangles)                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                    %
% EToV: triangle (nelem,3);                                 %
% VX(:,1): longitude (nnode,1);                             %
% VX(:,2): laritude (nnode,1);                              %
% B: depth (nnode,1);                                       %
% opedat.nbdv: Open boundary nodes ID                       %
% outfiname: mesh name                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali August 2018 ali.abdolali@noaa.gov       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('write WW3 msh') ;

node(:,1)=(1:length(VX(:,1)));
node(:,2)=VX(:,1);
node(:,3)=VX(:,2);
node(:,4)=B;

fid = fopen(outfiname,'w');
outname = outfiname ;
disp( title )  ;
fprintf(fid,'%s\n', '$MeshFormat');
fprintf(fid,'%s\n', '2 0 8');
fprintf(fid,'%s\n', '$EndMeshFormat');
fprintf(fid,'%s\n', '$Nodes');
fprintf(fid,'%d\n', length(node(:,1)));
   for i=1:length(node(:,1))
        fprintf(fid,['%d %s %5.5f %s %5.5f %s %5.5f\n'], node(i,1),'', VX(i,1),'',VX(i,2),'',B(i,1));
    end
fprintf(fid,'%s\n', '$EndNodes');
fprintf(fid,'%s\n', '$Elements');
fprintf(fid,'%d\n', length(EToV(:,1))+opedat.neta);
   m=0;
     for i=1:opedat.neta
       m=m+1;
       fprintf(fid,['%d %s %d %s %d %s %d %s %d %s %d\n'], m,'',15,'',2,'',0,'',0,'',opedat.nbdv(i));
     end

    for i=1:length(EToV(:,1))
       m=m+1;
      fprintf(fid,['%d %s %d %s %d %s %d %s %d %s %d %s %d %s %d %s %d\n'], m,'',2,'',3,'',0,'',i,'',0,'',EToV(i,1),'',EToV(i,2),'',EToV(i,3));
    end
fprintf(fid,'%s', '$EndElements');
fclose(fid) ;
return
