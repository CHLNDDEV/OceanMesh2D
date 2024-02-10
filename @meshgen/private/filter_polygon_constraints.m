function [pfix, egfix] = filter_polygon_constraints(pfix, egfix, ibouboxes, box_number)
%FILTER_CONSTRAINTS Removes edge constraints based on bounding box criteria.
%   This function filters out edge constraints (egfix) where at least one
%   endpoint (from pfix) is outside the specified bounding box (ibouboxes)
%   for the given box_number. It also adjusts the edges based on nested
%   bounding boxes, if applicable.
% remove bars if one point is outside
node1=pfix(egfix(:,1),:);
node2=pfix(egfix(:,2),:);
iboubox = ibouboxes{box_number};
% to enlarge or shrink the box, you must make sure bbox is equi
% spaced
[ty,tx]=my_interpm(iboubox(:,2),iboubox(:,1),100/111e3);
iboubox = [tx,ty]; 
buffer_size = 1.0;
iboubox(:,1) = buffer_size*iboubox(:,1)+(1-buffer_size)*mean(iboubox(1:end-1,1));
iboubox(:,2) = buffer_size*iboubox(:,2)+(1-buffer_size)*mean(iboubox(1:end-1,2));
inside_node1 = inpoly(node1,iboubox(1:end-1,:)) ;
inside_node2 = inpoly(node2,iboubox(1:end-1,:)) ;
inside = inside_node1 .* inside_node2;
% Get all points inside inner boxes and consider these outside for
% all the nested boxes.
for bn = box_number+1:length(ibouboxes)
    % enlarge nest
    iboubox = ibouboxes{bn}(1:end-1,:);
    [ty,tx]=my_interpm(iboubox(:,2),iboubox(:,1),100/111e3);
    iboubox = [tx,ty];
    iboubox(:,1) = 1.25*iboubox(:,1)+(1-1.25)*mean(iboubox(1:end-1,1));
    iboubox(:,2) = 1.25*iboubox(:,2)+(1-1.25)*mean(iboubox(1:end-1,2));
    inside_node1 = inpoly(node1,iboubox);
    inside_node2 = inpoly(node2,iboubox);
    inside2 = inside_node1 .* inside_node2;
    inside(find(inside2)) = false;
end
egfix(~inside,:) = [];
tegfix=egfix';
uid = unique(tegfix(:));
tuid = length(uid);
% remove nfix outside iboubox
if tuid > 0
    % remove pfix if outside domain
    pfix = pfix(uid,:);
end
egfix = renumberEdges(egfix);

end
