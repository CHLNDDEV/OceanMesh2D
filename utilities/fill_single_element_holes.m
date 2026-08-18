function [p_out,t_out,added_tris] = fill_single_element_holes(p,t)
%FILL_SINGLE_ELEMENT_HOLES Detect and fill single-triangle holes in a mesh.
%
% Syntax
%   [p_out,t_out,added] = fill_single_element_holes(p,t)
%
% Inputs
%   p  - (N x 2) node coordinates [lon, lat] or projected XY
%   t  - (M x 3) element connectivity (1-based indices into p)
%
% Outputs
%   p_out     - updated node coordinates (same as input unless fixmesh renumbers)
%   t_out     - updated element connectivity with hole-filling triangles added
%   added_tris- number of triangles added
%
% Notes
% - Uses existing OceanMesh2D utilities (extdom_edges2, extdom_polygon, fixmesh)
% - A "single element hole" is detected as an internal boundary loop with exactly
%   three unique vertices (i.e., a triangular void). The triangle formed by those
%   three nodes is added if it does not already exist.
%
% Dependencies
%   utilities/extdom_edges2.m
%   utilities/extdom_polygon.m
%   utilities/fixmesh.m
%
% Keith Roberts et al. OceanMesh2D (GPLv3)

% Defensive checks
if isempty(p) || isempty(t)
    p_out = p; t_out = t; added_tris = 0; return;
end

% Compute boundary edges (outer boundary + any holes)
[bnde,~] = extdom_edges2(t,p);
if isempty(bnde)
    p_out = p; t_out = t; added_tris = 0; return;
end

% Organize boundary edges into polygon loops (counter-clockwise ordering)
% poly:      cell array of (K_i x 2) coordinates following each boundary loop
% poly_idx:  cell array of (K_i x 1) point indices corresponding to poly
[poly,poly_idx] = extdom_polygon(bnde,p,0);

% Prepare set of existing triangles (sorted rows for membership tests)
if isempty(t)
    ts = zeros(0,3);
else
    ts = sort(t,2);
end

add_list = zeros(0,3);

for i = 1:numel(poly_idx)
    idx = poly_idx{i};
    if isempty(idx)
        continue
    end
    % Remove NaN separators if present
    idx = idx(~isnan(idx));
    if numel(idx) < 3
        continue
    end
    % Unique vertices in first-seen order (boundary walk order)
    u = unique(idx,'stable');
    if numel(u) ~= 3
        % Only fill triangular holes
        continue
    end
    tri_seq = u(:)'; % 1x3 in boundary order

    % Ensure positive (counter-clockwise) area for consistency
    A = poly_area(p(tri_seq,1), p(tri_seq,2));
    if A < 0
        tri_seq = fliplr(tri_seq);
    end

    % Add only if not already present
    tri_sorted = sort(tri_seq);
    if isempty(ts)
        exists = false;
    else
        exists = ismember(tri_sorted, ts, 'rows');
    end
    if ~exists
        add_list(end+1,:) = tri_seq; %#ok<AGROW>
        ts(end+1,:) = tri_sorted;    %#ok<AGROW>
    end
end

if ~isempty(add_list)
    t_new = [t; add_list];
    % Clean and relabel if needed
    [p_out,t_out] = fixmesh(p,t_new);
    added_tris = size(add_list,1);
else
    p_out = p; t_out = t; added_tris = 0;
end

end

function area = poly_area(x,y)
% Shoelace formula for polygon area (triangle variant)
% Assumes x,y are column vectors with 3 entries
area = 0.5 * ( x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)) );
end
