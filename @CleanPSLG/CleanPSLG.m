classdef CleanPSLG
    properties
        Vertices
        Segments
        Tolerance = 0.001;
        AngleThreshold = 10;
    end
    
    methods
        function obj = CleanPSLG(vertices, segments, tol, angle_thresh)
            obj.Vertices = vertices;
            obj.Segments = segments;
            if nargin > 2, obj.Tolerance = tol; end
            if nargin > 3, obj.AngleThreshold = angle_thresh; end
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
    end
end
