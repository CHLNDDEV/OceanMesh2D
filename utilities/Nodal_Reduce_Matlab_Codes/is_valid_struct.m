% IS_VALID_STRUCT - Determine if the input structure is "valid"
%
%     IS_VALID_STRUCT determines whether or not the input
%     structure is valid or invalid according to the 
%     description of FEM_GRID_STRUCT.  A valid FEM_GRID_STRUCT
%     contains (atleast) the following (NOT EMPTY) fields:
%
%       .name - the domain name of the FEM grid
%       .e    - node connectivity list (linear, triangular) [ne x 3 double]
%       .x    - x-horizontal node coordinates               [nn x 1 double]
%       .y    - y-horizontal node coordinates               [nn x 1 double]
%       .z    - bathymetry list                             [nn x 1 double]
%       .bnd  - boundary segment list                       [nb x 2 double]
%
%     As a (not very rigorous) check, the maximum node number in the 
%     element and boundary lists must not exceed the length of the node
%     list, and the node list ane bathymetry lists must be the same length.
%
%     There is no way to guarantee that the data in the fields actually
%     all came from the same domain (mesh).  However, if the 
%     FEM_GRID_STRUCT was output from LOADGRID or LOADG, then it is guaranteed.
%
% CALL: errflag=is_valid_struct(fem_grid_struct)
%       
% Written by : Brian O. Blanton
% Summer 1997
%
function errflag=is_valid_struct(fem_grid_struct)

errflag=0;

% Make sure input argument is actually a structure
%
if ~isstruct(fem_grid_struct)
   disp('    Argument to IS_VALID_STRUCT must be a structure.');return
end

% now, make sure the structure contains the minimum FEM_GRID_STRUCT fields
%
if ~isfield(fem_grid_struct,'name')
   disp('    Domain name field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'e')
   disp('    Element list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'x')
   disp('    X-node list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'y')
   disp('    Y-node list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'z')
   disp('    Bathymetry list field not part of fem_grid_struct');return
elseif ~isfield(fem_grid_struct,'bnd')
   disp('    Boundary list field not part of fem_grid_struct');return
end

% now, make sure the minimum FEM_GRID_STRUCT fields are NOT EMPTY
%
if isempty(fem_grid_struct.name)
   disp('    Domain name field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.e)
   disp('    Element list field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.x)
   disp('    X-node list field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.y)
   disp('    Y-node list field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.z)
   disp('    Bathymetry list field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.bnd)
   disp('    Boundary list field in fem_grid_struct is EMPTY');return
end

% Check the max node numbers in .bnd and .e fields
%
len_x=length(fem_grid_struct.x);
len_y=length(fem_grid_struct.y);
len_z=length(fem_grid_struct.z);
max_e=max(max(fem_grid_struct.e));
max_b=max(max(fem_grid_struct.bnd));

if(len_x~=len_y)
   disp('    Lengths of .x and .y fields NOT equal');return
elseif(len_x~=len_z)
   disp('    Lengths of .x, .y and .z fields NOT equal');return
elseif(max_e>len_x)
   disp('    Maximum node number in .e field exceeds length of .xy field');return
elseif(max_b>len_x)
   disp('    Maximum node number in .bnd field exceeds length of .xy field');return
end

errflag=1;

% These are non-fatal warnings!
% 
%if(max_e<len_xy)
%   disp('OPNML: non-fatal warning!!')
%   disp('    Maximum node number in .e field LESS THAN length of .xy field');return
%elseif(max_b<len_xy)
%   disp('OPNML: non-fatal warning!!')
%   disp('    Maximum node number in .bnd field LESS THAN length of .xy field');return
%end

%
%        Brian O. Blanton
%        Curr. in Marine Science
%        15-1A Venable Hall
%        CB# 3300
%        Uni. of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@marine.unc.edu
%
%        Summer 1997
