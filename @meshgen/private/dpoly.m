function d = dpoly(p,obj,box_vec,project)
% INPUTS:
%
% p: are the points to calculate the signed distance for.
%
% obj: is the meshgen object input from which we use...
%      outer is the bounding polygon (as a NaN-delimited vector).
%      bbox is the bounding box as a NaN-delimited vector (in CCW order).
%      inpoly_flip is whether to flip (1) the inpoly result or not (0).
%
% box_vec: is an integer (vector) denoting which box the calculation is in.
%          the calculation is looped over each box
% project: indicates to use the m_map toolbox to calculate the distance in
%            a projected space
%
% OUTPUTS:
%
% d is the signed distance from point p to closest point on polygon outer.
% (d is  negative if inside the bounded polygon outer and positive if outside it).
%
%     by Keith Roberts and William Pringle 2017-2018.
%        updated for clarity by Keith Roberts June 2019.
%%
% Generally this is always set in recent versions of the codebase.
if nargin < 4
    project = 1;
end

% If box_vec is not passed, construct it. 
if nargin == 2
    box_vec = 1:length(obj.bbox);
    % If box_vec is not passed and a projection switch is passed (for
    % reprojection/deluany_elim of points only)
elseif isempty(box_vec)
    box_vec = 1:length(obj.bbox);
end

% Initialize d with some positive number larger than geps
d = 0*p(:,1) + 1;

% For each box
for box_num = box_vec
    %  object ISN'T a cell
    if ~iscell(obj.outer)
        outer = obj.outer;
        inpoly_flip = obj.inpoly_flip;
        inside = true(size(d));
    else
        % geodata object IS a cell
        outer = obj.outer{box_num};
        inpoly_flip  = obj.inpoly_flip{box_num};
        % if it's not the outermost box
        if box_num > 1
            pt=obj.boubox{box_num};
            ee=Get_poly_edges(pt);
            inside = inpoly(p,pt,ee) ;
            clearvars pt ee
        else
            % It is the outermost box so all points are inside
            % this is enforced by construction of the boxes!
            inside = true(size(d));
        end
        % Get all points inside inner boxes and consider these outside for
        % all the nested boxes.
        for bn = box_num+1:length(obj.bbox)
            pt=obj.boubox{bn};
            ee=Get_poly_edges(pt);
            inside2 = inpoly(p,pt,ee) ;
            inside(inside2) = false;
        end
    end
    % for the points inside the box, calculate the nearest distance to
    % outer
    if sum(inside)~=0
        [~,d_l] = WrapperForKsearch(obj.anno{box_num},obj.annData{box_num},...
            p(inside,:),project);
    end
    
    %% Doing the signed calculation
    % The meshing domain is defined as the intersection of two areas:
    %
    % in_boubox are the points inside the bounding box 
    % in_outer are the points inside the outer meshing domain.
    %
    % "in" is defined as the intersection of the in = in_boubox && in_outer
    %
    % Note: the boubox is pre-pended when forming outer.
    firstNaN = find(isnan(outer(:,1)),1,'first') ;
    in_boubox = inpoly(p(inside,:),outer(1:firstNaN-1,:)) ;
    
    % Check if the points are inside the meshing domain (omega)
    edges = Get_poly_edges( outer );
    if sum(inside)~=0
        in_outer    = inpoly(p(inside,:),outer,edges);
        if inpoly_flip
            in_outer = ~in_outer;
        end
    else
        in_outer    = d_l*0;
    end
    
    % d is signed negative if inside and vice versa.
    d_l = (-1).^( in_outer & in_boubox).*d_l;
    
    if sum(inside)==0; return; end
    
    d(inside) = d_l;
end
end
