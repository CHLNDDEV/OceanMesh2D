function edges = Get_line_edges(nodes)
%GET_LINE_EDGES Generate edges from a sequence of nodes in a closed loop.
%   edges = GET_LINE_EDGES(nodes) takes a list of nodes (vertices) as input
%   and generates a list of edges where each node is connected to the next,
%   forming a closed loop. The output 'edges' is an Mx2 matrix, where M is
%   the number of edges, and each row represents an edge defined by the
%   indices of its two endpoints.

    % Number of nodes
    nNodes = length(nodes);

    % Generate edges connecting each node to the next
    edges = [(1:nNodes-1)', (2:nNodes)'];

    % Add the edge connecting the last node back to the first
    edges(end+1, :) = [nNodes, 1];
end
