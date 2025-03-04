function [pfix, egfix] = filter_polygon_constraints(pfix, egfix, ibouboxes, box_number)
%FILTER_POLYGON_CONSTRAINTS Removes edges and points outside a buffered bounding polygon.

    % Suppress polyshape inconsistency warning
    warning('off', 'MATLAB:polyshape:repairedBySimplify');

    % Extract bounding polygon and apply buffer
    iboubox = ibouboxes{box_number};
    buffered_poly = polybuffer(polyshape(iboubox), 100 / 111000); % ~100m buffer in degrees
    [tx, ty] = boundary(buffered_poly);
    iboubox = [tx, ty];

    % Extract node positions from edges
    node1 = pfix(egfix(:,1), :);
    node2 = pfix(egfix(:,2), :);

    % Identify edges where at least one endpoint is inside
    inside = inpolygon(node1(:,1), node1(:,2), iboubox(:,1), iboubox(:,2)) | ...
             inpolygon(node2(:,1), node2(:,2), iboubox(:,1), iboubox(:,2));

    % Process nested bounding boxes (if applicable)
    for bn = box_number+1:length(ibouboxes)
        nested_iboubox = ibouboxes{bn};
        nested_buffered_poly = polybuffer(polyshape(nested_iboubox), 125 / 111000); % ~125m buffer
        [tx, ty] = boundary(nested_buffered_poly);
        nested_iboubox = [tx, ty];

        % Identify edges where both endpoints are inside the nested box (remove these)
        inside_nested = inpolygon(node1(:,1), node1(:,2), nested_iboubox(:,1), nested_iboubox(:,2)) & ...
                        inpolygon(node2(:,1), node2(:,2), nested_iboubox(:,1), nested_iboubox(:,2));
        inside(inside_nested) = false;
    end

    % Retain only edges that pass filtering
    egfix = egfix(inside, :);

    % Retain only necessary points
    if ~isempty(egfix)
        uid = unique(egfix(:)); % Unique point indices used in edges
        pfix = pfix(uid, :); % Retain only those points
        egfix = renumberEdges(egfix); % Renumber edges to reflect reduced point set
    end

    % Restore warnings
    warning('on', 'MATLAB:polyshape:repairedBySimplify');

end