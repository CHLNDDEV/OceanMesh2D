function fem = laplace_smooth(fem_struct,nodenums,iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAPLACE_SMOOTH performs a Laplacian smoothing, aka "springing", on the 
%  mesh. This is achieved by moving each non-boundary point to the average
%  location of the x,y (,z) coordinates of the nodes connected to the given
%  point. The process is repeated according to the 'iter' input parameter.
%  The smoothing may also be performed on a limited portion of the mesh 
%  using the 'nodenums' input parameter.
%
% Usage -- fem = laplace_smooth(fem_struct,nodenums,iter);
%                ** nodenums & iter are optional **
%
% Variables
%  fem = updated finite element mesh
%  fem_struct = original finite element mesh
%  nodenums = node numbers of the nodes to be smoothed
%              (default is all nodes)
%  iter = number of iterations of smoothing to perform
%              (default is 5 iterations)
%
% Filename: laplace_smooth.m
% Created by: Ben Holladay
% Date: May 8, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parse input parameters. Throw exceptions for mal-formed input.
%
% Desired input variables: fem_struct, (nodenums), (iter)
%
% Defaults and valid flags
iterDefault = 5;
% nodenumsDefault = all non-boundary nodes
%
% Switch on number of input parameters.
switch nargin 
   case 0
      error('Not enough input parameters: laplace_smooth() must have a fem_struct.');
   case 1
      iter = iterDefault;
      nodenums = 1:length(fem_struct.x);      
   case 2
      iter = iterDefault;
      test = (floor(nodenums) == nodenums) & (nodenums > 0) & ...
         (nodenums <= length(fem_struct.x));
      if ~all(test)
         error('Input parameter ''nodenums'' must be a vector of valid node numbers.');
      end
   case 3
      iterValid = (floor(iter) == iter) && (iter > 0);
      if ~((floor(iter) == iter) && (iter > 0))
         error('Input parameter ''iter'' must be a positive integer');
      end
      test = (floor(nodenums) == nodenums) & (nodenums > 0) & ...
         (nodenums <= length(fem_struct.x));
      if ~all(test)
         error('Input parameter ''nodenums'' must be a vector of valid node numbers.');
      end
   otherwise
      error('Too many input parameters.');
end

% Set flag if their is no z variable.
try 
   z = fem_struct.z;
   zflag = true;
catch
   zflag = false;
end


% Initialize basic fem_struct variables.
enodes = fem_struct.e;
x = fem_struct.x;
y = fem_struct.y;
bnd = unique(fem_struct.bnd);
np = size(x,1);

% Nodes to be moved.
moveable = setdiff(1:length(x),bnd);
move = intersect(moveable,nodenums);
conflict = setdiff(nodenums,moveable);
if ~isempty(conflict)
   display('Boundary nodes will not be moved.');
end

% Spring the number of times given by iter.
for k = 1:iter
   % Create sparse matrix of neighbor points x,y,&z coordinates.
   tempx = sparse(enodes,enodes(:,[2,3,1]),x(enodes),np,np);
   tempy = sparse(enodes,enodes(:,[2,3,1]),y(enodes),np,np);
   if zflag
      tempz = sparse(enodes,enodes(:,[2,3,1]),z(enodes),np,np);
   end

   % Calculate the number of neighbors.
   numnghb = sparse(enodes,enodes(:,[2,3,1]),1,np,np);
   numnghb = sum(numnghb);

   % Find the new location.
   tempx = sum(tempx) ./ numnghb;
   tempy = sum(tempy) ./ numnghb;

   % Move nodes.
   x(move) = tempx(move);
   y(move) = tempy(move);

   % Find new location move z coordinate.
   if zflag
      tempz = sum(tempz) ./ numnghb;
      z(move) = tempz(move);
   end
end

% Create output fem_struct.
fem = fem_struct;
fem.x = x;
fem.y = y;
if zflag
   fem.z = z;
end
if is_valid_struct(fem_struct)
    fem = el_areas(fem);
end

return