function fem=springs(fem,iter);
% SPRING is legacy code that moves internal points of a finite 
%  element mesh such that the mesh quality is improved.  The movement
%  of nodes is accomplished with a Laplacian smoothing action, so 
%  that each point not located along the boundary is moved toward
%  the center of mass of the polygon formed by the adjacent
%  triangles. The process is repeated according to the setting 
%  of the ITER variable.
%
%       Calls laplace_smooth.m with specified input.
% FileName:  springs.m
% Written By:  unknown
% Date Last Modified:
%    May 9, 2007 -- Chris Massey, NRL Code 7322, Stennis Spc Cnt., MS
%                   Deleted code and replaced with a call to 
%                   laplace_smooth.m.
%    May 15, 2007 -- Chris Massey, -- Fixed bug, should send over the
%                   node list not the element list.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
if nargin==1,
    iter=5;
end

fem = laplace_smooth(fem,1:length(fem.x),iter);

return
