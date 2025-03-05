classdef CleanPSLG
    properties
        Vertices
        Segments
        Tolerance = 0.001;
        AngleThreshold = 10;
        DistanceThreshold = 0.001;  % Distance threshold for nearby edges
        PreservedSegments = [];     % Logical array marking segments to preserve
    end

    methods
        function obj = CleanPSLG(vertices, segments, tol, angle_thresh, dist_thresh)
            obj.Vertices = vertices;
            obj.Segments = segments;
            if nargin > 2, obj.Tolerance = tol; end
            if nargin > 3, obj.AngleThreshold = angle_thresh; end
            if nargin > 4, obj.DistanceThreshold = dist_thresh;
            else obj.DistanceThreshold = obj.Tolerance; end

            % Initialize preserved segments array
            obj.PreservedSegments = false(size(segments, 1), 1);
        end

        function obj = mergeVertices(obj)
            [obj.Vertices, obj.Segments] = obj.fix_geo(obj.Vertices, obj.Segments, obj.Tolerance);
        end

        function obj = dropIntersectingEdges(obj)
            obj.Segments = obj.drop_intersecting_edges(obj.Vertices, obj.Segments);
        end

        function obj = pruneEncroachingEdges(obj)
            obj.Segments = obj.drop_encroaching_edges(obj.Vertices, obj.Segments, obj.AngleThreshold, obj.Tolerance);
        end

        function obj = dropNearbyEdges(obj)
            obj.Segments = obj.drop_nearby_edges(obj.Vertices, obj.Segments, obj.DistanceThreshold, obj.PreservedSegments);
        end

        function obj = preserveLongChains(obj, minChainLength)
            % Preserve contiguous chains of segments longer than minChainLength
            % These chains will be protected during dropNearbyEdges

            % Identify chains of connected segments
            chains = obj.identifyContiguousChains();

            % Calculate the total length of each chain and mark for preservation
            numChains = length(chains);
            chainLengths = zeros(numChains, 1);
            preservedChains = false(numChains, 1);

            for i = 1:numChains
                chainLengths(i) = obj.calculateChainLength(chains{i});
                % Convert length to kilometers (assuming coordinates are in degrees)
                chainLengthKm = chainLengths(i) * 111; % Approximate conversion from degrees to km
                if chainLengthKm >= minChainLength
                    preservedChains(i) = true;
                end
            end

            % Mark segments in long chains as preserved
            obj.PreservedSegments = false(size(obj.Segments, 1), 1);
            for i = 1:numChains
                if preservedChains(i)
                    obj.PreservedSegments(chains{i}) = true;
                end
            end

            fprintf('Identified %d chains, preserving %d chains longer than %.1f km\n', ...
                numChains, sum(preservedChains), minChainLength);
            fprintf('Preserving %d segments out of %d total segments\n', ...
                sum(obj.PreservedSegments), length(obj.PreservedSegments));
        end

        function obj = deleteShortChains(obj, maxChainLength)
            % Delete chains of segments shorter than maxChainLength (in kilometers)
            % Identify chains of connected segments
            chains = obj.identifyContiguousChains();

            % Calculate the total length of each chain and identify short ones to delete
            numChains = length(chains);
            chainLengths = zeros(numChains, 1);
            shortChains = false(numChains, 1);
            segmentsToDelete = [];

            for i = 1:numChains
                chainSegs = chains{i};
                chainLengths(i) = obj.calculateChainLength(chainSegs);
                % Convert length to kilometers (assuming coordinates are in degrees)
                chainLengthKm = chainLengths(i) * 111; % Approximate conversion from degrees to km

                % Mark short chains for deletion
                if chainLengthKm < maxChainLength
                    shortChains(i) = true;
                    segmentsToDelete = [segmentsToDelete; chainSegs];
                end
            end

            % Remove the short chain segments
            if ~isempty(segmentsToDelete)
                keepMask = true(size(obj.Segments, 1), 1);
                keepMask(segmentsToDelete) = false;
                obj.Segments = obj.Segments(keepMask, :);

                % Also update the preserved segments array if it exists
                if ~isempty(obj.PreservedSegments)
                    obj.PreservedSegments = obj.PreservedSegments(keepMask);
                end
            end

            fprintf('Identified %d chains, deleted %d chains shorter than %.1f km\n', ...
                numChains, sum(shortChains), maxChainLength);
            fprintf('Deleted %d segments out of %d original segments (%.1f%%)\n', ...
                length(segmentsToDelete), length(segmentsToDelete) + size(obj.Segments, 1), ...
                100 * length(segmentsToDelete) / (length(segmentsToDelete) + size(obj.Segments, 1)));
        end

        function chains = identifyContiguousChains(obj)
            % Identify all contiguous chains of segments
            segments = obj.Segments;
            vertices = obj.Vertices;

            % Create vertex adjacency lists
            n = size(vertices, 1);
            vertexToSegments = cell(n, 1);

            for i = 1:size(segments, 1)
                v1 = segments(i, 1);
                v2 = segments(i, 2);
                vertexToSegments{v1} = [vertexToSegments{v1}, i];
                vertexToSegments{v2} = [vertexToSegments{v2}, i];
            end

            % Find connected components using depth-first search
            visited = false(size(segments, 1), 1);
            chains = {};

            for i = 1:size(segments, 1)
                if ~visited(i)
                    % Start a new chain
                    chain = [];
                    stack = i;

                    % DFS to find all connected segments
                    while ~isempty(stack)
                        current = stack(end);
                        stack(end) = [];

                        if ~visited(current)
                            visited(current) = true;
                            chain = [chain; current];

                            % Find connected segments through shared vertices
                            v1 = segments(current, 1);
                            v2 = segments(current, 2);

                            % Add connected segments to stack
                            for connectedSegIdx = [vertexToSegments{v1}, vertexToSegments{v2}]
                                if ~visited(connectedSegIdx)
                                    stack = [stack, connectedSegIdx];
                                end
                            end
                        end
                    end

                    chains{end+1} = chain;
                end
            end
        end

        function length = calculateChainLength(obj, chainIndices)
            % Calculate the total length of a chain of segments
            length = 0;
            for i = 1:numel(chainIndices)
                segIdx = chainIndices(i);
                v1 = obj.Vertices(obj.Segments(segIdx, 1), :);
                v2 = obj.Vertices(obj.Segments(segIdx, 2), :);
                length = length + norm(v2 - v1);
            end
        end

        function plotPSLG(obj)
            figure;
            subplot(1,2,1);
            scatter(obj.Vertices(:,1), obj.Vertices(:,2), 10, 'b', 'filled'); hold on;
            for i = 1:size(obj.Segments,1)
                pts = obj.Vertices(obj.Segments(i,:),:);
                plot(pts(:,1), pts(:,2), 'k-');
            end
            title('Processed PSLG'); xlabel('Longitude'); ylabel('Latitude'); axis equal;
        end

        function obj = cleanUnusedVertices(obj)
            % Remove vertices not referenced by any segments and renumber segments accordingly
            % Find all vertex indices that are used in segments
            usedIndices = unique(obj.Segments(:));
            numOriginalVertices = size(obj.Vertices, 1);

            % Create a mapping from old indices to new indices
            maxIndex = max(max(obj.Segments));
            vertexMap = zeros(maxIndex, 1);
            vertexMap(usedIndices) = 1:length(usedIndices);

            % Update the segments with new vertex indices
            obj.Segments = vertexMap(obj.Segments);

            % Keep only used vertices
            obj.Vertices = obj.Vertices(usedIndices, :);

            fprintf('Removed %d unused vertices (%.1f%% reduction)\n', ...
                numOriginalVertices - length(usedIndices), ...
                100 * (1 - length(usedIndices)/numOriginalVertices));
        end

    end

    methods (Static)
        function [newVertices, newSegments] = fix_geo(vertices, segments, tol)
            % Merge close vertices
            n = size(vertices,1);
            parent = (1:n)';
            function r = find_parent(i)
                while parent(i) ~= i, parent(i) = parent(parent(i)); i = parent(i); end
                r = i;
            end
            function union(i, j)
                ri = find_parent(i);
                rj = find_parent(j);
                if ri ~= rj, parent(rj) = ri; end
            end
            for i = 1:n-1
                for j = i+1:n
                    if norm(vertices(i,:) - vertices(j,:)) < tol
                        union(i,j);
                    end
                end
            end
            mapping = arrayfun(@find_parent, (1:n)');
            [uniqueMappings, ~, newMapping] = unique(mapping);
            newVertices = vertices(uniqueMappings, :);
            newSegments = arrayfun(@(x) newMapping(x), segments);
            newSegments = unique(sort(newSegments, 2), 'rows');
        end

        function newSegments = drop_intersecting_edges(vertices, segments)
            % Drop intersecting edges
            nseg = size(segments,1);
            drop = false(nseg,1);
            for i = 1:nseg
                for j = i+1:nseg
                    if any(segments(i,:) == segments(j,:))
                        continue;
                    end
                    p = vertices(segments(i,1),:);
                    r = vertices(segments(i,2),:);
                    q = vertices(segments(j,1),:);
                    s = vertices(segments(j,2),:);
                    if CleanPSLG.segments_intersect(p, r, q, s)
                        drop(i) = true; drop(j) = true;
                    end
                end
            end
            newSegments = segments(~drop,:);
        end

        function flag = segments_intersect(p, r, q, s)
            o1 = CleanPSLG.orientation(p, r, q);
            o2 = CleanPSLG.orientation(p, r, s);
            o3 = CleanPSLG.orientation(q, s, p);
            o4 = CleanPSLG.orientation(q, s, r);
            flag = (o1 ~= o2) && (o3 ~= o4);
        end

        function o = orientation(p, q, r)
            val = (q(2)-p(2))*(r(1)-q(1)) - (q(1)-p(1))*(r(2)-q(2));
            o = (val > 0) - (val < 0);
        end

        function newSegments = drop_encroaching_edges(vertices, segments, angle_thresh_deg, tol)
            nseg = size(segments,1);
            drop = false(nseg,1);
            cos_thresh = cosd(angle_thresh_deg);
            for i = 1:nseg
                if drop(i), continue; end
                for j = i+1:nseg
                    if drop(j), continue; end
                    if any(segments(i,:) == segments(j,:))
                        continue;
                    end
                    v1 = vertices(segments(i,:),:);
                    v2 = vertices(segments(j,:),:);
                    unit1 = (v1(2,:) - v1(1,:)) / norm(v1(2,:) - v1(1,:));
                    unit2 = (v2(2,:) - v2(1,:)) / norm(v2(2,:) - v2(1,:));
                    if abs(dot(unit1, unit2)) > cos_thresh
                        drop(i) = true;
                    end
                end
            end
            newSegments = segments(~drop,:);
        end

        function newSegments = drop_nearby_edges(vertices, segments, dist_thresh, preservedSegments)
            % Drop edges that are close to other edges - optimized version
            % Only considers perpendicular distances between edges
            % Respects segments marked for preservation
            nseg = size(segments, 1);
            if nseg < 2
                newSegments = segments;
                return;
            end

            % If preservedSegments not provided, assume no segments are preserved
            if nargin < 4 || isempty(preservedSegments)
                preservedSegments = false(nseg, 1);
            end

            % Calculate segment midpoints and bounding boxes for spatial binning
            midpoints = zeros(nseg, 2);
            bbox = zeros(nseg, 4); % [minX, minY, maxX, maxY]

            for i = 1:nseg
                v1 = vertices(segments(i,1), :);
                v2 = vertices(segments(i,2), :);
                midpoints(i, :) = (v1 + v2) / 2;
                bbox(i, :) = [min(v1(1), v2(1)), min(v1(2), v2(2)), ...
                    max(v1(1), v2(1)), max(v1(2), v2(2))];
            end

            % Create spatial grid for efficient lookup
            % Determine grid dimensions
            domain = [min(bbox(:,1)), min(bbox(:,2)), max(bbox(:,3)), max(bbox(:,4))];
            gridSize = max(1, ceil(dist_thresh * 10)); % Adjust based on dist_thresh

            numCellsX = max(1, ceil((domain(3) - domain(1)) / gridSize));
            numCellsY = max(1, ceil((domain(4) - domain(2)) / gridSize));

            % Assign segments to grid cells
            cellAssignment = cell(numCellsX, numCellsY);
            for i = 1:nseg
                % Determine which cells this segment overlaps
                minCellX = max(1, min(numCellsX, floor((bbox(i,1) - domain(1)) / gridSize) + 1));
                minCellY = max(1, min(numCellsY, floor((bbox(i,2) - domain(2)) / gridSize) + 1));
                maxCellX = max(1, min(numCellsX, ceil((bbox(i,3) - domain(1)) / gridSize)));
                maxCellY = max(1, min(numCellsY, ceil((bbox(i,4) - domain(2)) / gridSize)));

                % Add segment to all overlapping cells
                for cx = minCellX:maxCellX
                    for cy = minCellY:maxCellY
                        cellAssignment{cx, cy} = [cellAssignment{cx, cy}; i];
                    end
                end
            end

            % Process segments
            drop = false(nseg, 1);

            % Fast check for potential nearby segments using grid
            for cx = 1:numCellsX
                for cy = 1:numCellsY
                    cellSegments = cellAssignment{cx, cy};
                    numCellSegs = length(cellSegments);

                    if numCellSegs < 2
                        continue;  % Skip cells with 0 or 1 segment
                    end

                    % Check all pairs of segments in this cell
                    for ii = 1:numCellSegs-1
                        i = cellSegments(ii);
                        if drop(i) || preservedSegments(i)
                            continue;  % Skip already dropped or preserved segments
                        end

                        for jj = ii+1:numCellSegs
                            j = cellSegments(jj);
                            if drop(j) || preservedSegments(j) || any(segments(i,:) == segments(j,:))
                                continue;  % Skip already dropped, preserved, or connected segments
                            end

                            % Quick bounding box check
                            if bbox(i,1) > bbox(j,3) + dist_thresh || ...
                                    bbox(j,1) > bbox(i,3) + dist_thresh || ...
                                    bbox(i,2) > bbox(j,4) + dist_thresh || ...
                                    bbox(j,2) > bbox(i,4) + dist_thresh
                                continue;  % Bounding boxes too far apart
                            end

                            % Detailed perpendicular distance check
                            v1 = vertices(segments(i,:), :);
                            v2 = vertices(segments(j,:), :);

                            % Get perpendicular distances
                            [perp_dist, is_perp] = CleanPSLG.perpendicular_segment_distance(v1, v2);

                            % Only consider perpendicular distances
                            if is_perp && perp_dist < dist_thresh
                                % Drop the longer segment (unless preserved)
                                len1 = norm(v1(2,:) - v1(1,:));
                                len2 = norm(v2(2,:) - v2(1,:));
                                if len1 > len2
                                    drop(i) = true;
                                    break;  % Break inner loop, segment i is dropped
                                else
                                    drop(j) = true;
                                end
                            end
                        end
                    end
                end
            end

            newSegments = segments(~drop,:);
        end

        function [dist, is_perpendicular] = perpendicular_segment_distance(v1, v2)
            % Calculate perpendicular distance between two line segments
            % Returns the minimum perpendicular distance and whether it actually occurs
            % v1 and v2 are 2x2 matrices containing start and end points of segments

            % Get segment vectors
            p1 = v1(1,:);
            p2 = v1(2,:);
            p3 = v2(1,:);
            p4 = v2(2,:);

            % Direction vectors of the segments
            dir1 = p2 - p1;
            dir2 = p4 - p3;

            % Initialize values
            dist = Inf;
            is_perpendicular = false;

            % Normalize directions
            len1 = norm(dir1);
            len2 = norm(dir2);

            if len1 < 1e-10 || len2 < 1e-10
                % One of the segments is almost a point
                return;
            end

            unit_dir1 = dir1 / len1;
            unit_dir2 = dir2 / len2;

            % Check if segments are nearly parallel (not perpendicular)
            dot_product = abs(dot(unit_dir1, unit_dir2));
            if dot_product > 0.9  % Approximately parallel if dot product > 0.9 (< 26 degrees difference)
                % Segments are nearly parallel - not considering them perpendicular
                return;
            end

            % For the first segment, calculate perpendicular distance to the second segment
            d1 = CleanPSLG.point_to_segment_perp_distance(p1, p3, p4);
            d2 = CleanPSLG.point_to_segment_perp_distance(p2, p3, p4);

            % For the second segment, calculate perpendicular distance to the first segment
            d3 = CleanPSLG.point_to_segment_perp_distance(p3, p1, p2);
            d4 = CleanPSLG.point_to_segment_perp_distance(p4, p1, p2);

            % Get the minimum perpendicular distance
            min_perp_dist = min([d1, d2, d3, d4]);

            % Only consider it a perpendicular case if the minimum distance is finite
            if isfinite(min_perp_dist)
                dist = min_perp_dist;
                is_perpendicular = true;
            end
        end

        function dist = point_to_segment_perp_distance(p, s1, s2)
            % Calculate perpendicular distance from point p to segment (s1,s2)
            % Returns Inf if projection is outside the segment

            v = s2 - s1;  % Segment vector
            w = p - s1;   % Vector from segment start to point

            % Check if we can project onto the segment
            c1 = dot(w, v);
            if c1 <= 0
                % p projects outside the segment, beyond s1
                dist = Inf;
                return;
            end

            c2 = dot(v, v);
            if c2 <= c1
                % p projects outside the segment, beyond s2
                dist = Inf;
                return;
            end

            % Calculate projection parameter and perpendicular distance
            b = c1 / c2;  % Projection parameter (0 to 1)
            pb = s1 + b * v;  % Projected point on segment
            dist = norm(p - pb);
        end

        function d = fast_segment_distance(v1, v2)
            % Optimized calculation of minimum distance between two line segments

            % First do quick endpoint checks
            d_endpoints = min([
                norm(v1(1,:) - v2(1,:)),
                norm(v1(1,:) - v2(2,:)),
                norm(v1(2,:) - v2(1,:)),
                norm(v1(2,:) - v2(2,:))
                ]);

            % Early exit if endpoints are close enough
            if d_endpoints < 1e-10
                d = 0;
                return;
            end

            % Get line segment vectors and cross-segment vector
            p1 = v1(1,:);
            p2 = v1(2,:);
            p3 = v2(1,:);
            p4 = v2(2,:);

            v13 = p1 - p3;
            v43 = p4 - p3;
            v21 = p2 - p1;

            % Check if either segment is actually a point
            d43 = sum(v43.^2);
            d21 = sum(v21.^2);

            % Handle special cases
            if d43 < 1e-10 && d21 < 1e-10  % Both segments are points
                d = norm(p1 - p3);
                return;
            end

            if d21 < 1e-10  % First segment is a point
                d = CleanPSLG.point_line_distance(p1, p3, p4);
                return;
            end

            if d43 < 1e-10  % Second segment is a point
                d = CleanPSLG.point_line_distance(p3, p1, p2);
                return;
            end

            % Compute the parameters of the closest points
            d1343 = sum(v13 .* v43);
            d4321 = sum(v43 .* v21);
            d1321 = sum(v13 .* v21);
            d4343 = d43;
            d2121 = d21;

            denom = d2121 * d4343 - d4321 * d4321;

            if abs(denom) < 1e-10  % Lines are parallel
                % Use perpendicular distance between parallel lines
                d_p1l2 = CleanPSLG.point_line_distance(p1, p3, p4);
                d = d_p1l2;
                return;
            end

            numer = d1343 * d4321 - d1321 * d4343;

            mua = numer / denom;
            mub = (d1343 + d4321 * mua) / d4343;

            % Clamp parameters to segment bounds
            mua = max(0, min(1, mua));
            mub = max(0, min(1, mub));

            % Compute the closest points on each segment
            pa = p1 + mua * v21;
            pb = p3 + mub * v43;

            % Return the distance between these points
            d = norm(pa - pb);
        end

        function d = segment_distance(v1, v2)
            % Calculate minimum distance between two line segments
            % v1 and v2 are 2x2 matrices containing start and end points of segments

            % Vector directions of the segments
            d1 = v1(2,:) - v1(1,:);
            d2 = v2(2,:) - v2(1,:);

            % Point-to-point distances
            d_p1p1 = norm(v1(1,:) - v2(1,:));
            d_p1p2 = norm(v1(1,:) - v2(2,:));
            d_p2p1 = norm(v1(2,:) - v2(1,:));
            d_p2p2 = norm(v1(2,:) - v2(2,:));

            % Point-to-line distances
            d_p1l2 = CleanPSLG.point_line_distance(v1(1,:), v2(1,:), v2(2,:));
            d_p2l2 = CleanPSLG.point_line_distance(v1(2,:), v2(1,:), v2(2,:));
            d_p3l1 = CleanPSLG.point_line_distance(v2(1,:), v1(1,:), v1(2,:));
            d_p4l1 = CleanPSLG.point_line_distance(v2(2,:), v1(1,:), v1(2,:));

            % Minimum distance is the smallest of all distances
            d = min([d_p1p1, d_p1p2, d_p2p1, d_p2p2, d_p1l2, d_p2l2, d_p3l1, d_p4l1]);
        end

        function d = point_line_distance(p, l1, l2)
            % Calculate perpendicular distance from point p to line through l1 and l2
            % Check if projection falls within line segment
            v = l2 - l1;
            w = p - l1;

            c1 = dot(w, v);
            if c1 <= 0
                d = norm(p - l1);
                return;
            end

            c2 = dot(v, v);
            if c2 <= c1
                d = norm(p - l2);
                return;
            end

            b = c1 / c2;
            pb = l1 + b * v;
            d = norm(p - pb);
        end


    end
end