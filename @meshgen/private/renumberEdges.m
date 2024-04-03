function edges = renumberEdges(edges)
%RENUMBEREDGES Renumber edges to maintain consistency after point modifications.
%   edges = RENUMBEREDGES(edges) updates the indices of edge endpoints in
%   the 'edges' matrix to reflect a new, sequential numbering based on
%   their occurrence in 'edges'. This is useful after points have been
%   filtered or reordered, ensuring that edge references are consistent
%   with the updated point list.

% Find unique indices of edge endpoints and create a direct mapping
idx = unique(edges(:));
for i = 1 : length(idx)
    map(idx(i),1) = i;
end

% renumber bnde to correspond with pfix
for i = 1 : size(edges,1)
    edges(i,1) = map(edges(i,1));
    edges(i,2) = map(edges(i,2));
end


end
