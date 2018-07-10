% EL_AREAS - compute triangular finite element areas
%
%         EL_AREAS(FEM_GRID_STRUCT) computes the areas for the 
%         elements of the FEM domain described in the 
%         structure FEM_GRID_STRUCT.  The function must
%         return a new structure, which is identical to the 
%         input structure with the element areas attached.
%         The element areas are contained in the field .ar,
%         so that the new structure now includes:
%
%         .ar   - element areas                    [ne x 1 double]  
%         
%
%         EL_AREAS optionally returns an index of element
%         numbers whose areas are negative (if any).
%         Negative element areas indicate clockwise elemental
%         node numbering, instead of the conventional counter-
%         clockwise numbering.
%
%  INPUT : fem_grid_struct - (from LOADGRID, see FEM_GRID_STRUCT)       
%           
% OUTPUT : new_struct (REQ) - new structure with areas
%          ineg       (OPT) - index of negative area elements
%
%   CALL : >>new_struct=el_areas(fem_grid_struct);
%    or 
%          >>[new_struct,ineg]=el_areas(fem_grid_struct);
%
% Written by : Brian O. Blanton 
% Summer 1997
%

function [ret_struct,ineg]=el_areas(fem_grid_struct)

% VERIFY INCOMING STRUCTURE
%
if ~is_valid_struct(fem_grid_struct)
   error('    Argument to EL_AREAS must be a valid fem_grid_struct.')
end

% NEED ONE or TWO OUT ARGS
%
if nargout~=1 & nargout~=2
   error('   EL_AREAS must have 1 or 2 output arguments.')
end

% BREAK DOWN INCOMING STRUCTURE
%
e=fem_grid_struct.e;
x=fem_grid_struct.x;
y=fem_grid_struct.y;

% COMPUTE GLOBAL DY
%
dy=[y(e(:,2))-y(e(:,3)) y(e(:,3))-y(e(:,1)) y(e(:,1))-y(e(:,2))];

% COMPUTE ELEMENTAL AREAS
%
AR=(x(e(:,1)).*dy(:,1)+x(e(:,2)).*dy(:,2)+x(e(:,3)).*dy(:,3))/2.;

%Create return structure and attach element areas to ret_struct
%
ret_struct=fem_grid_struct;
ret_struct.ar=AR;

if nargout==2
   % ANY NEGATIVE OR ZERO AREAS ?
   %
   ineg=find(AR<=0);
end

%
%        Brian O. Blanton
%        Curriculum in Marine Sciences
%        Ocean Processes Numerical Modeling Laboratory
%        15-1A Venable Hall
%        CB# 3300
%        Uni. of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@marine.unc.edu
%
%        October 1995
%
