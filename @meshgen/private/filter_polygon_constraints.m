function [pfix, egfix] = filter_polygon_constraints(pfix, egfix, ibouboxes, box_number)
%FILTER_POLYGON_CONSTRAINTS Filters edge constraints based on bounding box criteria.
% Uses MATLAB's polybuffer to properly buffer the bounding box.
% Skips the inside check if the original and buffered polygons have similar areas.
    % Suppress polyshape inconsistency warning
    warning('off', 'MATLAB:polyshape:repairedBySimplify');

    % Extract node positions for each edge
    node1 = pfix(egfix(:,1), :);
    node2 = pfix(egfix(:,2), :);

    % Convert bounding box to polyshape
    iboubox = ibouboxes{box_number};
    iboubox_shape = polyshape(iboubox(:,1), iboubox(:,2));

    % Define buffer size in degrees (~100m, adjust if necessary)
    buffer_distance = 100 / 111000; % Convert meters to degrees

    % Apply buffer using polybuffer
    buffered_iboubox = polybuffer(iboubox_shape, buffer_distance);

    % Compute areas of original and buffered bounding boxes
    original_area = area(iboubox_shape);
    buffered_area = area(buffered_iboubox);

    % Define area similarity threshold (5% tolerance)
    area_tolerance = 0.05;
    
    % If areas are similar, skip the inside check
    if abs(buffered_area - original_area) / original_area < area_tolerance
        return; % Keep all edges without filtering
    end

    % Extract updated boundary
    [tx, ty] = boundary(buffered_iboubox);
    iboubox = [tx, ty];

    % Check if nodes are inside the buffered bounding box
    inside_node1 = inpolygon(node1(:,1), node1(:,2), iboubox(:,1), iboubox(:,2));
    inside_node2 = inpolygon(node2(:,1), node2(:,2), iboubox(:,1), iboubox(:,2));

    % Keep edges where at least one endpoint is inside
    inside = inside_node1 | inside_node2;

    % Handle nested bounding boxes with same logic
    for bn = box_number+1:length(ibouboxes)
        nested_iboubox = ibouboxes{bn};
        nested_iboubox_shape = polyshape(nested_iboubox(:,1), nested_iboubox(:,2));

        % Apply a larger buffer for nested boxes
        nested_buffered_iboubox = polybuffer(nested_iboubox_shape, buffer_distance * 1.25);

        % Compute areas for nested box
        nested_area = area(nested_iboubox_shape);
        nested_buffered_area = area(nested_buffered_iboubox);

        % Skip inside check if areas are similar
        if abs(nested_buffered_area - nested_area) / nested_area < area_tolerance
            continue;
        end

        % Extract boundary
        [tx, ty] = boundary(nested_buffered_iboubox);
        nested_iboubox = [tx, ty];

        % Identify edges where both endpoints are inside the nested box
        inside_node1 = inpolygon(node1(:,1), node1(:,2), nested_iboubox(:,1), nested_iboubox(:,2));
        inside_node2 = inpolygon(node2(:,1), node2(:,2), nested_iboubox(:,1), nested_iboubox(:,2));

        % Remove edges fully inside the nested box
        inside_nested = inside_node1 & inside_node2;
        inside(inside_nested) = false;
    end

    % Remove edges that are fully outside all bounding boxes
    egfix = egfix(inside, :);

    % Ensure pfix only contains necessary points
    if ~isempty(egfix)
        tegfix = egfix';
        uid = unique(tegfix(:));
        pfix = pfix(uid, :);
        egfix = renumberEdges(egfix);
    end

    warning('on', 'MATLAB:polyshape:repairedBySimplify');

end
