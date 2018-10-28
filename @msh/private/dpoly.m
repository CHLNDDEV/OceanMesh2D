function d = dpoly(p,obj,box_vec)
% d = dpoly(p,obj,box_vec)
%
% p are the mesh points
% and obj is the meshgen object input where we use...
% outer is the bounding polygon
% bbox is the bounding box
% inpoly_flip is whether to flip the inpoly result or not

% d is the distance from point p to closest point on polygon outer.
% (d is  negative if inside the bounded polygon outer and positive if outside it).
% by Keith Roberts and William Pringle 2017-2018.

%% Doing the distance calc
if nargin == 2
    box_vec = 1:length(obj.bbox);
elseif isempty(box_vec)
    box_vec = 1:length(obj.bbox);
end
% initialize d with some positive number larger than geps
d = ones(length(p),1);
for box_num = box_vec
    if ~iscell(obj.outer)
        outer = obj.outer;
        pv1 = outer;
        pv1(isnan(obj.outer(:,1)),:) = [];
        inpoly_flip = obj.inpoly_flip;
        bbox = obj.bbox;
        %h0_l = obj.h0;
        inside = true(size(d));
    else
        outer = obj.outer{box_num};
        pv1 = outer;
        pv1(isnan(obj.outer{box_num}(:,1)),:) = [];
        inpoly_flip  = obj.inpoly_flip{box_num};
        bbox = obj.bbox{box_num};
        if box_num > 1
            inside = (p(:,1) >= bbox(1,1) & p(:,1) <= bbox(1,2) & ...
                p(:,2) >= bbox(2,1) & p(:,2) <= bbox(2,2) );
        else
            inside = true(size(d));
        end
        for bn = box_num+1:length(obj.bbox)
            % Get all points inside inner boxes
            inside2 = (p(:,1) >= obj.bbox{bn}(1,1) & ...
                p(:,1) <= obj.bbox{bn}(1,2) & ...
                p(:,2) >= obj.bbox{bn}(2,1) & ...
                p(:,2) <= obj.bbox{bn}(2,2) );
            inside(inside2) = false;
        end
    end
    
    if sum(inside)~=0
        [~,d_l] = WrapperForKsearch(pv1', p(inside,:)');
    end
    
    %% Doing the inpoly check
    % If inpoly m file we need to get the edges to pass to it to avoid the
    % issues of NaNs
    edges = Get_poly_edges( outer );
    if sum(inside)~=0
        in    = inpoly(p(inside,:),outer,edges);
    else
        in    = d_l*0; 
    end
    
    % d is negative if inside polygon and vice versa.
    if inpoly_flip
        d_l = (-1).^(~in).*d_l;
    else
        d_l = (-1).^( in).*d_l;
    end
    
    if sum(inside)==0; return; end;
        
    % IF OUTSIDE BUT APPEARS INSIDE
    bad = find((p(inside,1) < bbox(1,1) | p(inside,1) > bbox(1,2) | ...
        p(inside,2) < bbox(2,1) | p(inside,2) > bbox(2,2)) & d_l < 0);
    if nnz(bad) > 0
        d_l(bad) = -d_l(bad);
    end
    d(inside) = d_l;
end
%EOF
end
