function [points_new, triangles_new] = collapseMinAngleTriangles(points, triangles, angleThreshold, eligibleMask)
    % Collapses triangles with minimum angle below threshold
    %
    % Inputs:
    %   points - Nx2 array of (x,y) coordinates
    %   triangles - Mx3 array of triangle connectivity (1-based indexing)
    %   angleThreshold - Minimum angle threshold in degrees (default: 6°)
    %   eligibleMask - [optional] Boolean mask of triangles eligible for collapsing
    %
    % Outputs:
    %   points_new - Updated points after collapsing
    %   triangles_new - Updated triangulation after collapsing
    
    if nargin < 3
        % Default: 6 degrees = 0.1 * 60 degrees (10% of equilateral)
        angleThreshold = 6.0;
    end
    
    % Calculate quality based on minimum angle
    [quality, minAngles] = minAngleMeshQuality(points, triangles);
    
    % Convert angle threshold to normalized quality
    qualityThreshold = angleThreshold / 60.0;
    
    % Find triangles below quality threshold
    poorQualityMask = quality < qualityThreshold;
    
    % Apply eligibility mask if provided
    if nargin > 3
        poorQualityMask = poorQualityMask & eligibleMask;
    end
    
    % Find poor quality triangle indices
    poorTriangles = find(poorQualityMask);
    
    % Copy inputs to outputs for modification
    points_new = points;
    triangles_new = triangles;
    
    % Process each poor quality triangle
    for i = 1:length(poorTriangles)
        triIdx = poorTriangles(i);
        
        % Skip if this triangle no longer exists
        if triIdx > size(triangles_new, 1)
            continue;
        end
        
        % Get vertices of the triangle
        v1 = triangles_new(triIdx, 1);
        v2 = triangles_new(triIdx, 2);
        v3 = triangles_new(triIdx, 3);
        
        % Calculate edge lengths
        edge1 = norm(points_new(v2,:) - points_new(v3,:)); % Edge opposite to v1
        edge2 = norm(points_new(v1,:) - points_new(v3,:)); % Edge opposite to v2
        edge3 = norm(points_new(v1,:) - points_new(v2,:)); % Edge opposite to v3
        
        % Find the shortest edge
        [~, minEdgeIdx] = min([edge1, edge2, edge3]);
        
        % Identify vertices to merge based on shortest edge
        switch minEdgeIdx
            case 1 % Merge v2 and v3
                vKeep = v2;
                vRemove = v3;
            case 2 % Merge v1 and v3
                vKeep = v1;
                vRemove = v3;
            case 3 % Merge v1 and v2
                vKeep = v1;
                vRemove = v2;
        end
        
        % Calculate midpoint
        midpoint = (points_new(vKeep,:) + points_new(vRemove,:)) / 2;
        
        % Update position of kept vertex to midpoint
        points_new(vKeep,:) = midpoint;
        
        % Update connectivity: replace all instances of vRemove with vKeep
        triangles_new(triangles_new == vRemove) = vKeep;
        
        % Remove degenerate triangles (those with repeated vertices)
        degenerate = false(size(triangles_new,1), 1);
        for j = 1:size(triangles_new,1)
            % Direct comparison is faster than unique for 3 elements
            v = triangles_new(j,:);
            if (v(1) == v(2)) || (v(2) == v(3)) || (v(1) == v(3))
                degenerate(j) = true;
            end
        end
        triangles_new(degenerate,:) = [];
    end
    
    % Remove unused vertices and reindex
    [points_new, triangles_new] = cleanMesh(points_new, triangles_new);
end

function [points_clean, triangles_clean] = cleanMesh(points, triangles)
    % Removes unused vertices and reindexes the triangulation
    
    % Create boolean mask of used vertices (faster than unique)
    usedMask = false(size(points,1), 1);
    usedMask(triangles(:)) = true;  % Mark all used vertices
    
    % Get the indices of used vertices
    usedVertices = find(usedMask);
    
    % Create mapping from old to new indices
    vertexMap = zeros(size(points,1), 1);
    vertexMap(usedVertices) = 1:length(usedVertices);
    
    % Apply mapping to get new triangulation
    triangles_clean = vertexMap(triangles);
    
    % Extract only used vertices
    points_clean = points(usedVertices,:);
end