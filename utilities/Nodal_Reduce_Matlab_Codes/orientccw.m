function fem = orientccw(fem_struct)
% orientccw checks to see if the vertices in the element
% list (fem_struct.e) are arranged in a CCW direction.
% If any element's vertices are not in a CCW direction
% that element's vertices are rearranged.
%
% Variables:
% fem_struct == opnml fem_grid_struct
% fem == new fem_grid_struct
%
% Routines Called:
% el_areas.m
%
% CALL:  fem = orientccw(fem_struct);
%
% Filename:   orientccw.m
% Written By: Chris Massey
% Date:       Aug. 19, 2003
% Modified:   Sept. 3, 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keywords:  counter clockwise, negative area, rearrange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fem = fem_struct;

[fem,ineg]=el_areas(fem_struct);


if isempty(ineg) == 0
    temp = fem.e(ineg,2);
    fem.e(ineg,2) = fem.e(ineg,3);
    fem.e(ineg,3) = temp;
    clear temp
    fem.ar = abs(fem.ar);
end


return
