% DETBNDY compute a boundary segment list for a FEM domain
%
% DETBNDY bnd=detbndy(e,outfile);
%         This function computes a boundary for the FEM domain
%         described a file containing element connectivity list (e).
%         It uses sparse matrix techniques to determine the element
%         edges on the boundary of the FEM domain.
%
% Input:  ele -  element list; 3 (.tri) or 4 (.ele) columns wide
%         file - name of output file (.bnd extension)
%                (*optional*)
% Output: bnd -  a 2-column list of boundary-node numbers, returned
%                to the local workspace
%
%         The output boundary list are pairs of node numbers, not 
%         coordinates, describing the edges of elements on the 
%         exterior of the domain, including islands.  The segments 
%         are not connected.
%
%         Call as: bnd=detbndy(e,outfile);
%
% Written by : Brian O. Blanton at The University of North Carolina 
%              at Chapel Hill, Mar 1995.
%
function bnd=detbndy(in,outfile)

% DEFINE ERROR STRINGS
err1=['Max input arguments to DETBNDY. Type "help detbndy"'];
err2=['Element list passed to DETBNDY does not have 3 or 4 columns'];
err3=str2mat(['??? Error using ==> detbndy'],...
             ['DETBNDY must have one output argument.'],...
             ['Call as: >>  bnd=detbndy(in);'],...
             [' ']);
% check argument list

if nargin>2
   error(err1);
end

if nargout~=1
   disp(err3);
   return
end
 
 
% Check size of element list
[nelems,ncol]=size(in);
if ncol < 3 | ncol > 4
   error(err2);
   return
elseif ncol==4
   in=in(:,2:4);
end

% Form (i,j) connection list from .ele element list
%
i=[in(:,1);in(:,2);in(:,3)];
j=[in(:,2);in(:,3);in(:,1)];

% Form the sparse adjacency matrix and add transpose.
%
n = max(max(i),max(j));
A = sparse(i,j,-1,n,n);
A = A + A';

% Consider only the upper part of A, since A is symmetric
% 
A=A.*triu(A);

% The boundary segments are A's with value == 1
%
B=A==1;

% Extract the row,col from B for the boundary list.
%
[ib,jb,s]=find(B);
bnd=[ib(:),jb(:)];

% Output .bnd list

if nargin==2
   fp=fopen(outfile,'w');
   for i=1:length(bnd)
      fprintf(fp,'%d %d\n', bnd(i,1),bnd(i,2));  
   end
end
return;
%
%        Brian O. Blanton
%        Curriculum in Marine Sciences
%        15-1A Venable Hall
%        CB# 3300
%        University of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@marine.unc.edu
%

