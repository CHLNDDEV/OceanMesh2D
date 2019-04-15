function [vtoe,nne]=VertToEle(t)
np = max(t(:)); ne = size(t,1);
nne = zeros(np,1);
for ie = +1 : ne                                                            %--Go through once and find max nnz.
    %for iv = +1 : +3
    %kjr loop unrolled.
    nm1=t(ie,1);
    nm2=t(ie,2);
    nm3=t(ie,3);
    %vrtx = t(ie,iv);
    nne(nm1,1) = nne(nm1,1) + 1;
    nne(nm2,1) = nne(nm2,1) + 1;
    nne(nm3,1) = nne(nm3,1) + 1;
    %end
end
mnz  = max(nne);                                                           %--max number of non-zeros
vtoe = zeros(mnz,np);                                                      %--vertex to element connectivity
nne = zeros(np,1);                                                         %--number of neighboring elements
for ie = +1 : ne
    nm1=t(ie,1);
    nm2=t(ie,2);
    nm3=t(ie,3);
    
    %kjr loop unrolled.
    nne(nm1,1) = nne(nm1,1) +1;
    nne(nm2,1) = nne(nm2,1) +1;
    nne(nm3,1) = nne(nm3,1) +1;
    
    vtoe(nne(nm1,1),nm1) = ie; 
    vtoe(nne(nm2,1),nm2) = ie; 
    vtoe(nne(nm3,1),nm3) = ie; 

    %vtoe(nm1,nne(nm1,1)) = ie;
    %vtoe(nm2,nne(nm2,1)) = ie;
    %vtoe(nm3,nne(nm3,1)) = ie;
end
nne  = nne';
%vtoe = vtoe';

end