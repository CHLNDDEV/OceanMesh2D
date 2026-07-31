function [quality, minAngles] = minAngleMeshQuality(points, triangles)
    % Calculates mesh quality based solely on minimum angle criterion
    %
    % Inputs:
    %   points - Nx2 array of (x,y) coordinates
    %   triangles - Mx3 array of triangle connectivity (1-based indexing)
    %   angleThreshold - Minimum angle threshold in degrees (default: 6°)
    %
    % Outputs:
    %   quality - Normalized quality values (0 to 1, 1 is best)
    %   minAngles - Minimum angle in each triangle (in degrees)
    
    if nargin < 3
        % Default: 6 degrees = 0.1 * 60 degrees (10% of equilateral)
        angleThreshold = 6.0;
    end
    
    numTriangles = size(triangles, 1);
    quality = zeros(numTriangles, 1);
    minAngles = zeros(numTriangles, 1);
    
    for i = 1:numTriangles
        v1 = triangles(i,1);
        v2 = triangles(i,2);
        v3 = triangles(i,3);
        
        % Triangle vertices
        p1 = points(v1,:);
        p2 = points(v2,:);
        p3 = points(v3,:);
        
        % Edge vectors
        e1 = p2 - p3; % opposite to v1
        e2 = p3 - p1; % opposite to v2
        e3 = p1 - p2; % opposite to v3
        
        % Edge lengths
        len1 = norm(e1);
        len2 = norm(e2);
        len3 = norm(e3);
        
        % Calculate angles using law of cosines
        angles = zeros(1,3);
        
        % Ensure we don't get NaN from floating point issues
        cos1 = (len2^2 + len3^2 - len1^2) / (2 * len2 * len3);
        cos1 = min(max(cos1, -1), 1);
        
        cos2 = (len1^2 + len3^2 - len2^2) / (2 * len1 * len3);
        cos2 = min(max(cos2, -1), 1);
        
        cos3 = (len1^2 + len2^2 - len3^2) / (2 * len1 * len2);
        cos3 = min(max(cos3, -1), 1);
        
        angles(1) = acosd(cos1);
        angles(2) = acosd(cos2);
        angles(3) = acosd(cos3);
        
        minAngle = min(angles);
        minAngles(i) = minAngle;
        
        % Normalize by 60 degrees (equilateral triangle)
        % For minimum angles near 0, this will give quality values near 0
        quality(i) = minAngle / 60.0;
    end
end