function [obj,ind,in] = extract_subdomain(obj,bou,varargin)
% [obj,ind,in] = extract_subdomain(obj,bou,varargin)
%
% Inputs:
% bou: a bbox, i.e.: [lon min, lon_max;
%                     lat_min, lat_max], or
%      a NaN-delimited polygon of the domain to extract
%
% Varargins...
% keep_inverse: = 0 [default] to get the sub-domain inside the bou polygon
%               = 1 to get the sub-domain outside the bou polygon
% keep_numbering: = 0 [default] to renumber and "fix" the mesh
%               = 1 to keep the original triangulation numbering without mesh "fix"
% centroid: = 0 [default] inpolygon test is based on whether all vertices
%              of the element are inside (outside) the bou polygon
%           = 1 inpolygon test is based on whether the element centroid
%              is inside (outside) the bou polygon
% min_depth = topobathymetry value to keep elements inside bou so that they have a
%             sufficiently deep value at their centroid
% max_depth = topobathymetry value to keep elements inside bou so that they have a
%             sufficiently shallow value at their centroid.
% nscreen: = 1 [default] display the notice to screen
%          = 0 do not display the notice to screen
%
% Outputs:
% obj: the subset mesh obj (only p and t, properties untouched)
% ind: an array of indices that can be used to map the mesh properties to
%      the output mesh subset with a subsequent call to "map_mesh_properties".
%  in: the logical array of elements within the given boundary
%      used to perform the subsetting 
%
keep_inverse = 0 ;
keep_numbering = 0 ;
centroid = 0 ;
min_depth = -99999;
max_depth = +99999;
nscreen = 1;
% Otherwise, name value pairs specified.
% Parse other varargin
for kk = 1:2:length(varargin)
    if strcmp(varargin{kk},'keep_inverse')
        keep_inverse = varargin{kk+1};
    elseif strcmp(varargin{kk},'keep_numbering')
        keep_numbering = varargin{kk+1};
    elseif strcmp(varargin{kk},'centroid')
        centroid = varargin{kk+1};
    elseif strcmp(varargin{kk},'nscreen')
        nscreen = varargin{kk+1};
    elseif strcmp(varargin{kk},'min_depth')
        min_depth = varargin{kk+1};
    elseif strcmp(varargin{kk},'max_depth')
        max_depth = varargin{kk+1};
    end
end

p = obj.p; t = obj.t; b = obj.b;
if nargin == 1 || (nargin == 2 && isempty(bou))
    plot(p(:,1),p(:,2),'k.');
    h = impoly;
    bou  = h.getPosition;
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
if min_depth ~= -99999 || max_depth ~= +99999
     bem = max(b(t),[],2);   % only trim when all vertices
     selected = bem > min_depth & bem < max_depth;
     in = logical(in .* selected);
end
if keep_inverse == 0
    t(~in,:) = [];
else
    t(in,:) = [];
end

if keep_numbering
   ind = unique(t(:));
   obj.p = p(ind,:); obj.t = t;
else
   % Remove uncessary vertices and reorder triangulation
   [p1,t1,ind] = fixmesh(p,t);
   % Put back into the msh obj
   obj.p = p1; obj.t = t1;
end
%
if nscreen
    % Displaying notice for mapping mesh properties
    disp('NOTICE: Only p and t have been subset.')
    disp('  To map mesh properties to the subset output the ind array and call: ')
    disp('  mesh_obj = map_mesh_properties(mesh_obj,''ind'',ind)')
end
% EOF
end
