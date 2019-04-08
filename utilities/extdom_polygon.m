function [poly,poly_idx,max_index,max_size] = extdom_polygon(bnde,pts,order,line,min_size)
% DESCRIPTION: Given a set of boundary edges of a singly- or multi-
%              polygonal region, organize them in a winding order.
%
% INPUTS:
%         bnde: the indices of each boundary edge as a nbnde x 2 matrix
%         pts:  the x,y locations of all the points in the region
%               stored as an np x 2 matrix.
%         order:the order in which the traversal takes place
%               counter-clockwise (0) or clockwise (1) or add a negative
%               sign to append NaNs in each cell.
%          line:if desired output will be polylines
%       min_size: if the lenght of the polygon is less than min_size
%                  points, throw it out. zero by default. 
% OUTPUTS:
%          poly: the boundary of each enclosing polygon sorted in winding-order
%                poly is returned as a cell-array of length number of polys.
%      poly_idx: indices of the polygon coordinates in the same format as
%                poly
%     max_index: is the index into poly that is the largest
%      max_size: is the size of the largest poly
%
% Last Edited:
% kjr,UND,CHL,2017
%   kjr,UND,CHL -->revised for massive speed improvements March 2018.
%
%                                           TRAVERSAL METHOD
% Pick any unvisited edge segment [v_start,v_next] and add these vertices to the polygon loop.
% Find the unvisited edge segment [v_i,v_j] that has either v_i = v_next or v_j = v_next and add the other vertex (the one not equal to v_next) to the polygon loop.
% Reset v_next as this newly added vertex, mark the edge as visited and continue from 2.
% Traversal is done when we get back to v_start.
% NOTE: that the signed area will be positive if the vertices are
% oriented counterclockwise, and will be negative if it is oriented clockwise

% NOTE: By flipping the edges left-to-right, we can easily see the nodal
% connectivity of the triangulation. However, since we are effectively
% adding new "edges", we must also quickly locate and flag the "flipped" edge
% in addition to the current edge under consideration. This is accomplished with the invmap
% array which allows one, given a vertex gid, to quickly find it's linear
% index inside the array bnde((boundary edges). This localizes the search to find
% the connectivity to continue the "walk" on the boundary making the
% calculation massively more efficient than compared to searching the entire bnde
% array with edge under consideration.

% if storing lines
if(nargin < 4); line = 0; end
if(nargin < 5); min_size = 1; end

bnde = [bnde; fliplr(bnde)];
bnde = sortrows(bnde,1);
ned = length(bnde);
an  = sign(order);
active = true(ned,1);

p = 0 ;
% given a vertex with num gid, return where it exists last lid.
for lid = 1 : ned
    gid = bnde(lid,1);
    invmap(gid) = lid ;
end

while any(active)
    p = p + 1;
    
    rn = find(active,1);
    
    temp  = pts(bnde(rn,:)',:);
    temp2 = bnde(rn,:)';
    
    v_start= bnde(rn,1);
    v_next = bnde(rn,2);
    
    active(rn) = false;
    
    % flag flipped edge too
    flipped = fliplr(bnde(rn,:));
    idx = invmap(flipped(1));
    st   = max((idx - 1),1);
    ed   = min(idx,ned);
    rn = find(flipped(2)==bnde(st:ed,2),1);
    
    active(rn+st-1) = false;
    
    k = 2 ;
    while v_next~=v_start
        
        % form local set to search for continuation of boundary walk.
        idx  = invmap(v_next);
        st   = max((idx - 1),1);
        ed   = min(idx,ned);
        
        r = find(v_next==bnde(st:ed,1) & active(st:ed),1);
        tsel = bnde(st+r-1,:); tsel_inv = fliplr(tsel) ;
        sel=tsel(tsel~=v_next);
        
        if(line)
            if(isempty(sel))
                break
            end
        end
        % store points.
        k = k + 1;
        temp(k,:) = pts(sel,:);
        temp2(k,:)= sel;
        
        active(r+st-1) = false;
        
        % flag flipped edge too
        idx = invmap(tsel_inv(1));
        st   = max((idx - 1),1);
        ed   = min(idx,ned);
        r = find(tsel_inv(2)==bnde(st:ed,2),1);
        active(r+st-1) = false;
        
        v_next = sel;
        
    end
    
    if length(temp > min_size)
        poly{p}     = temp;
        poly_idx{p} = temp2;
    else
        continue 
    end
    if an < 0
        poly{p}(end+1,:) = [NaN NaN];
        poly_idx{p}(end+1,:) = NaN;
    end
    [area]=parea(poly{p}(:,1),poly{p}(:,2));
    if order==0 % ccw
        if sign(area)<0
            poly{p} = flipud(poly{p});
            poly_idx{p} = flipud(poly_idx{p});
        end
    else % cw
        if sign(area)>0
            poly{p} = flipud(poly{p});
            poly_idx{p} = flipud(poly_idx{p});
        end
    end
end
[max_size, max_index] = max(cellfun('size', poly, 1));
end
% helper function, computes area of polygon
function [area]=parea(x,y)
n    = length(x);
xp   = [x; x(1)];
yp   = [y; y(1)];
area = 0;
for i = 1:n
    area = area + det([xp(i), xp(i+1); yp(i), yp(i+1)]);
end
area = 1/2*area;
end