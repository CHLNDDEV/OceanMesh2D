function [edges]=renumberEdges(edges)

% create mapping to points
idx = unique(edges(:));
for i = 1 : length(idx)
    map(idx(i),1) = i;
end

% renumber bnde to correspond with pfix
for i = 1 : length(edges)
    edges(i,1) = map(edges(i,1));
    edges(i,2) = map(edges(i,2));
end


end