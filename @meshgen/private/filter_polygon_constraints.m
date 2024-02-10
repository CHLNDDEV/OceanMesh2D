function [pfix, egfix] = filter_constraints(pfix, egfix, ibouboxes, box_number)
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


% function [pfix, egfix] = filter_polygon_constraints(pfix, egfix, ibouboxes, box_number)
% %FILTER_CONSTRAINTS Removes edges based on bounding box constraints.
% %   [pfix, egfix] = FILTER_CONSTRAINTS(pfix, egfix, ibouboxes, box_number)
% %   adjusts the fixed points (pfix) and edge constraints (egfix) for a
% %   given bounding box (specified by box_number within ibouboxes). Edges
% %   are removed if either of their endpoints lies outside the main bounding
% %   box or inside any nested bounding boxes.
% 
% % Extract endpoints of each edge
% node1 = pfix(egfix(:,1), :);
% node2 = pfix(egfix(:,2), :);
% 
% % Main bounding box for the current box_number
% iboubox = ibouboxes{box_number};
% 
% % Interpolate to refine the bounding box resolution and adjust dimensions
% [ty, tx] = my_interpm(iboubox(:,2), iboubox(:,1), 100/111e3); % Custom interpolation
% iboubox = [tx, ty]; % Updated bounding box after interpolation
% 
% % Adjust bounding box by applying a buffer size to modify its scale
% buffer_size = 1.0; % No scaling applied if buffer_size is 1
% iboubox(:,1) = buffer_size * iboubox(:,1) + (1 - buffer_size) * mean(iboubox(1:end-1,1));
% iboubox(:,2) = buffer_size * iboubox(:,2) + (1 - buffer_size) * mean(iboubox(1:end-1,2));
% 
% % Determine if nodes are inside the main bounding box
% inside_node1 = inpoly(node1, iboubox(1:end-1, :));
% inside_node2 = inpoly(node2, iboubox(1:end-1, :));
% inside = inside_node1 & inside_node2; % Logical AND to find edges fully inside
% 
% % Exclude edges based on nested boxes
% for bn = box_number+1:length(ibouboxes)
%     % Process each nested bounding box
%     iboubox = ibouboxes{bn}(1:end-1,:);
%     [ty, tx] = my_interpm(iboubox(:,2), iboubox(:,1), 100/111e3);
%     iboubox = [tx, ty]; % Update bounding box after interpolation
% 
%     % Scale and adjust the nested bounding box
%     scale_factor = 1.25;
%     iboubox(:,1) = scale_factor * iboubox(:,1) + (1 - scale_factor) * mean(iboubox(1:end-1,1));
%     iboubox(:,2) = scale_factor * iboubox(:,2) + (1 - scale_factor) * mean(iboubox(1:end-1,2));
% 
%     % Check if nodes are inside this nested box
%     inside_node1 = inpoly(node1, iboubox);
%     inside_node2 = inpoly(node2, iboubox);
%     inside2 = inside_node1 & inside_node2; % Edges inside nested box
% 
%     % Mark edges inside nested boxes as outside
%     inside(inside2) = false;
% end
% 
% % Remove edges not meeting the inside criteria
% egfix(~inside, :) = [];
% 
% % Renumber edges to maintain consistency
% egfix = renumberEdges(egfix); % Assumes renumberEdges is a custom function
% 
% % Filter pfix based on unique edge endpoints
% uniqueEndpoints = unique(egfix(:));
% if ~isempty(uniqueEndpoints)
%     pfix = pfix(uniqueEndpoints, :);
% end
% 
% end
