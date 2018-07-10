% IS_VALID_STRUCT2 - Determine if the input structure is "valid"
%
%     IS_VALID_STRUCT2 determines whether or not the input
%     structure contains additional field data beyond the 
%     minimum valid fields checked for by  IS_VALID_STRUCT.
%     Recall that a valid FEM_GRID_STRUCT contains (atleast) 
%     the following (NOT EMPTY) fields:
%
%       .name,.e,.x,.y,.z,.bnd  
%
%     IS_VALID_STRUCT2 verifies the existence of .A, .B, .A0,
%     and .T which are filled by BELINT, and of .ar filled by
%     EL_AREAS.  These additional FEM fields are needed by 
%     OPNML/MATLAB element-finding, basis  and interpotion 
%     routines.
%
% CALL: errflag=is_valid_struct2(fem_grid_struct)
%       
% Written by : Brian O. Blanton
% Summer 1998
%
function errflag=is_valid_struct2(fem_grid_struct)

errflag=0;

% Make sure input argument is actually a structure
%
if ~isstruct(fem_grid_struct)
   disp('    Argument to IS_VALID_STRUCT2 must be a structure.');return
end

% now, make sure the structure contains the additional fields,
% as above
if ~isfield(fem_grid_struct,'A')
   disp('    "A" field not part of fem_grid_struct; run BELINT');return
elseif ~isfield(fem_grid_struct,'B')
   disp('    "B" field not part of fem_grid_struct; run BELINT');return
elseif ~isfield(fem_grid_struct,'A0')
   disp('    "A0" field not part of fem_grid_struct; run BELINT');return
elseif ~isfield(fem_grid_struct,'T')
   disp('    "T" field not part of fem_grid_struct; run BELINT');return
elseif ~isfield(fem_grid_struct,'ar')
   disp('    "ar" field not part of fem_grid_struct; run ');return
end

% now, make sure these additional fields are NOT EMPTY
%
if isempty(fem_grid_struct.A)
   disp('    "A" field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.B)
   disp('    "B" field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.A0)
   disp('    "A0" field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.T)
   disp('    "T" field in fem_grid_struct is EMPTY');return
elseif isempty(fem_grid_struct.ar)
   disp('    "ar" field in fem_grid_struct is EMPTY');return
end

errflag=1;

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
%        Summer 1998
