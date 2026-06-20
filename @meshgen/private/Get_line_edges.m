function edges = Get_line_edges(nodes)
%GET_LINE_EDGES Generate edges from a sequence of nodes.
%   edges = GET_LINE_EDGES(nodes) takes a list of nodes (vertices) as input
%   and generates a list of edges where each node is connected to the next.
%   If the nodes form a closed loop, the last node is connected back to the first.
%   The output 'edges' is an Mx2 matrix, where M is the number of edges,
%   and each row represents an edge defined by the indices of its two endpoints.
%
%   Assumes 'nodes' is an Nx2 array (or Nx3 for 3D), where each row is a coordinate.

    % Number of nodes
    nNodes = size(nodes, 1);

    % Generate edges connecting each node to the next
    edges = [(1:nNodes-1)', (2:nNodes)'];

    % Check if the shape is a closed loop
    if all(nodes(1, :) == nodes(end, :))
        % Add the edge connecting the last node back to the first
        edges(end+1, :) = [nNodes, 1];
    end
end
