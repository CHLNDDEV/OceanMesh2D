function [obj,ind] = extract_subdomain(obj,bou,keep_inverse,centroid,threshold,nscreen)
% [obj,ind] = extract_subdomain(obj,bou,keep_inverse,centroid,threshold,nscreen)
% 
% Inputs:
% bou: a bbox, i.e.: [lon min, lon_max;
%                     lat_min, lat_max], or    
%      a NaN-delimited polygon of the domain to extract
% keep_inverse: = 0 [default] to get the sub-domain inside the bou polygon
%               = 1 to get the sub-domain outside the bou polygon
% centroid: = 0 [default] inpolygon test is based on whether all vertices
%              of the element are inside (outside) the bou polygon
%           = 1 inpolygon test is based on whether the element centroid
%              is inside (outside) the bou polygon
% threshold = bathymetry value to keep elements inside bou so that they have a
%             sufficiently deep "threshold" at they're centroid
% nscreen: = 1 [default] display the notice to screen
%          = 0 do not display the notice to screen
% 
% Outputs:
% obj: the subset mesh obj (only p and t, properties untouched)
% ind: an array of indices that can be used to map the mesh properties to
% the output mesh subset with a subsequent call to "map_mesh_properties". 
%
p = obj.p; t = obj.t; b = obj.b;
if nargin == 1 || (nargin == 3 && isempty(bou))
    plot(p(:,1),p(:,2),'k.');
    h = impoly;
    bou  = h.getPosition;
end
if nargin < 3 || isempty(keep_inverse)
    keep_inverse = 0;
end
if nargin < 4 || isempty(centroid)
    centroid = 0;
end
if nargin < 5 || isempty(threshold)
    threshold = -9999;
end
if nargin < 6 || isempty(nscreen)
    nscreen = 1;
end
% converting bbox to polygon
if size(bou,1) == 2
     bou = [bou(1,1) bou(2,1);
            bou(1,1) bou(2,2); 
            bou(1,2) bou(2,2);
            bou(1,2) bou(2,1); 
            bou(1,1) bou(2,1)];
end
nans = false;
% check whether there are nans in the bou
if sum(isnan(bou(:,1))) > 0; nans = true; end
if nans
    edges = Get_poly_edges(bou);
end
if centroid
    bxyc = baryc(obj);
    if ~nans
        in = inpoly(bxyc,bou);
    else
        in = inpoly(bxyc,bou,edges);
    end
else
    bxy1 = p(t(:,1),:); bxy2 = p(t(:,2),:); bxy3 = p(t(:,3),:); 
    if ~nans
        in1 = inpoly(bxy1,bou); 
        in2 = inpoly(bxy2,bou); 
        in3 = inpoly(bxy3,bou);
    else
        in1 = inpoly(bxy1,bou,edges); 
        in2 = inpoly(bxy2,bou,edges); 
        in3 = inpoly(bxy3,bou,edges);
    end
    in = in1 & in2 & in3;
end
if threshold ~= -9999
     bem = max(b(t),[],2);   % only trim when all vertices
     overland = bem < -threshold; % of element are above the threshold.
     in = logical(in .* overland);
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
    disp('  mesh_obj = map_mesh_properties(mesh_obj,ind)')
end
% EOF
end
