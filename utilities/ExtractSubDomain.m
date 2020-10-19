function [obj,ind] = ExtractSubDomain(obj,bou,keep_inverse,centroid,nscreen)
% [obj,ind] = ExtractSubDomain(obj,bou,keep_inverse,centroid,nscreen)
% 
% Inputs:
% bou: a polygon or a bbox to extract sub-domain (no NaNs allowed)
% keep_inverse: = 0 [default] to get the sub-domain inside the bou polygon
%               = 1 to get the sub-domain outside the bou polygon
% centroid: = 0 [default] inpolygon test is based on whether all vertices
%              of the element are inside (outside) the bou polygon
%           = 1 inpolygon test is based on whether the element centroid
%              is inside (outside) the bou polygon
% nscreen: = 1 [default] display the notice to screen
%          = 0 do not display the notice to screen
% 
% Outputs:
% obj: the subset mesh obj (only p and t, properties untouched)
% ind: an array of indices that can be used to map the mesh properties to
% the output mesh subset with a subsequent call to "mapMeshProperties". 
%
p = obj.p; t = obj.t; 
if nargin == 1 || (nargin == 3 && isempty(bou))
    plot(p(:,1),p(:,2),'k.');
    h = impoly;
    bou  = h.getPosition;
end
if nargin < 3
    keep_inverse = 0;
end
if nargin < 4
    centroid = 0;
end
if nargin < 5
    nscreen = 1;
end
if size(bou,1) == 2
     bou = [bou(1,1) bou(2,1);
            bou(1,1) bou(2,2); ...
            bou(1,2) bou(2,2);
            bou(1,2) bou(2,1); ...
            bou(1,1) bou(2,1)];
end
if centroid
    bxyc = baryc(obj);
    in = inpoly(bxyc,bou);
else
    bxy1 = p(t(:,1),:); bxy2 = p(t(:,2),:); bxy3 = p(t(:,3),:); 
    in1 = inpoly(bxy1,bou); in2 = inpoly(bxy2,bou); in3 = inpoly(bxy3,bou);
    in = in1 & in2 & in3;
end
if keep_inverse == 0
    t(~in,:) = [];
else
    t(in,:) = [];
end
% Remove uncessary vertices and reorder triangulation
[p1,t,ind] = fixmesh(p,t);
% Put back into the msh obj
obj.p = p1; obj.t = t;  
%
if nscreen
    % Displaying notice for mapping mesh properties
    disp('NOTICE: Only p and t have been subset.')
    disp('  To map mesh properties to the subset output the ind array and call: ')
    disp('  mesh_obj = mapMeshProperties(mesh_obj,ind)')
end
% EOF
end