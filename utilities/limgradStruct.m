function [ffun,flag] = limgradStruct(ny,xeglen,yeglen,ffun,fdfdx,imax)
%LIMGRAD impose "gradient-limits" on a function defined over
%an undirected graph.
%   [FNEW] = LIMGRADStruct(NY,EGLEN,FFUN,DFDX,IMAX) computes a
%   "gradient-limited" function FNEW on the undirected structured graph
%   that has Ny rows when supplied with a fixed edgelength EGLEN. 
%   It uses an 8 node stencil to relax the gradients. 
%   FFUN must be numbered such that rows vary the fastest.
%   Gradients are limited over the graph edges, such that
%
%       ABS(FNEW(N2)-FNEW(N1)) <= ELEN(II) * DFDX,
%
%   where N1=EDGE(II,1) and N2=EDGE(II,2) are the two nodes
%   in the II-TH edge. An iterative algorithm is used, swee-
%   ping over an "active-set" of graph edges until converge-
%   nce is achieved. A maximum of IMAX iterations are done.
%
%   [FNEW,FLAG] = LIMGRADStruct(...) also returns a boolean FLAG,
%   with FLAG=TRUE denoting convergence.


%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 18/04/2017
% ---------------------
%         Modified for a structred graph with eight node stencil 
%         Last updated: 24/06/2017
%         Email       : krober10@nd.edu
%         Keith Roberts, 2017.
% ---------------------
%             Modified to have spatially variable fdfdx 
%             Last updated: 27/04/2019
%             Keith Roberts, 2019

if length(fdfdx)==1 
  fdfdx = ffun*0 + fdfdx ; 
end
%----------------------------- ASET=ITER if node is "active"
aset = zeros(size(ffun,1),1) ;

%----------------------------- exhaustive 'til all satisfied
ftol = min(ffun) * sqrt(eps) ;

% ------------------------ calculate hypotenuse eglen length
dgeglen = sqrt(xeglen.^2+repmat(yeglen,1,ny).^2); 
eglen   = 0.5*(xeglen+repmat(yeglen,1,ny)); 

rm = zeros(9,1); 
for iter = +1 : imax
    
    %------------------------- find "active" nodes this pass
    aidx = find(aset == iter - 1) ;
    
    if (isempty(aidx)), break; end
    %------------------------- reorder => better convergence
    [aval,idxx] = sort(ffun(aidx)) ;
    
    aidx = aidx(idxx);
    
    %------------------------- speed up a little by preallocating occasionally
    npos = zeros(9,1); elens = zeros(9,1); 
    
    %------------------------- visit adj. edges and set DFDX
    for i = 1 : length(aidx)
        % ----- map doubly index to singly indexed
        inod = aidx(i);
        ipos = 1 + floor((inod-1)/ny); 
        jpos = inod - (ipos - 1)*ny;
        
        % ------ use 8 edge stencil for higher order gradients
        nn=1;
        npos(nn) =  inod;                        nn=nn+1;
        npos(nn) =  ipos*ny    + jpos;           nn=nn+1;%--- nnod of right adj
        npos(nn) = (ipos-2)*ny + jpos;           nn=nn+1;%--- nnod of left adj
        npos(nn) = (ipos-1)*ny + min(jpos+1,ny); nn=nn+1;%--- nnod of above adj
        npos(nn) = (ipos-1)*ny + max(jpos-1,1);  nn=nn+1;%--- nnod of below adj
        npos(nn) = (ipos*ny)   +  max(jpos-1,1); nn=nn+1;%--- nnod of right bot diag adj
        npos(nn) =  (ipos*ny)  + min(jpos+1,ny); nn=nn+1;%--- nnod of right top diag adj
        npos(nn) = (ipos-2)*ny +  min(jpos+1,ny);nn=nn+1;%--- nnod of left top diag adj
        npos(nn) = (ipos-2)*ny + max(jpos-1,1);          %--- nnod of left bot diag adj
        
        
        % ------ populate elens
        nn = 1;
        elens(nn) = eglen(jpos);  nn=nn+1;
        elens(nn) = xeglen(jpos); nn=nn+1;
        elens(nn) = xeglen(jpos); nn=nn+1;
        elens(nn) = yeglen; nn=nn+1;
        elens(nn) = yeglen; nn=nn+1;
        elens(nn) = dgeglen(jpos); nn=nn+1;
        elens(nn) = dgeglen(jpos); nn=nn+1; 
        elens(nn) = dgeglen(jpos); nn=nn+1; 
        elens(nn) = dgeglen(jpos); 
        
        %----- handle boundary vertex adjs.
        rm = npos <= 0 | npos > size(ffun,1);
        npos(rm) = [];
        elens(rm)= [];
        
        for ne = 2 : length(npos)
            nod1 = npos(1);
            nod2 = npos(ne);
            elen = elens(ne);
            
            %----------------- calc. limits about min.-value
            if (ffun(nod1) > ffun(nod2))
                
                fun1 = ffun(nod2) ...
                    + elen * fdfdx(nod2) ;
                
                if (ffun(nod1) > fun1+ftol)
                    ffun(nod1) = fun1;
                    aset(nod1) = iter;
                end
                
            else
                
                fun2 = ffun(nod1) ...
                    + elen * fdfdx(nod2) ;
                
                if (ffun(nod2) > fun2+ftol)
                    ffun(nod2) = fun2;
                    aset(nod2) = iter;
                end
                
            end
        end
    end
    rm = rm*0; 
    flag = (iter < imax) ;
    
end



