function [etoe,idx] = EleToEle(t)
%--if an element shares an edge then it is a neighbor 
%--returns element-to-element connectivity in crs format. 
%--triangle ie is connected to etoe(idx(ie):idx(ie+1)-1,2) triangles
% kjr 2017, fixed August 2017 for self adj triangles (i.e., disjoint
% elements). 
t = sort(t,2);
edges = [t(:,[1 2]);t(:,[1 3]);t(:,[2 3])];
nt = size(t,1);
trinum = repmat((1:nt)',3,1);
% use sort by converting edges to bytes for speed.
[~,tags]=sort(edges*[2^31;1]);
edges=edges(tags,:); 
%[edges,tags] = sortrows(edges);
trinum = trinum(tags);
k = find(all(diff(edges,1)==0,2));
etoe=trinum([k,k+1]);
[~,dmy1]=sort(etoe(:,1));[~,dmy2]=sort(etoe(:,2)); 
% use sort by converting edges to bytes for speed
temp=[etoe(dmy1,:); fliplr(etoe(dmy2,:)); (1:nt)',(1:nt)'];
[~,tags]=sort(temp*[2^31;1]); 
etoe=temp(tags,:); 
%etoe=sortrows([etoe(dmy1,:); fliplr(etoe(dmy2,:)); (1:nt)',(1:nt)']);%added self adjs. 
idx = find(diff(etoe(:,1))==1); idx=[idx;length(etoe)]; idx = [1;idx+1];                                              
end