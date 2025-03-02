function polygon_struct = Read_shapefile(finputname, ~, bbox, h0, boubox, plot_on, ~)

%% Initialize Data Structures
SG = []; % Store read polygons
loop = 1; minus = 0;
tolerance = 1e-5; % **Define closure tolerance** (adjustable)

if bbox(1,2) > 180 && bbox(1,1) < 180, loop = 2; end
if all(bbox(1,:) > 180), minus = 1; end

%% Read Shapefile Data
for fname = finputname
    for nn = 1:loop
        bboxt = bbox';
        if loop == 2
            bboxt(2-((nn-1)*1),1) = 180 - (nn-1)*360;
        end
        if minus, bboxt(:,1) = bboxt(:,1) - 360; end

        try
            S = shaperead(fname{1}, 'BoundingBox', bboxt);
            S = rmfield(S, setdiff(fieldnames(S), {'X', 'Y'}));
            disp('Read shapefile with shaperead');
        catch
            S = m_shaperead(fname{1}, reshape(bboxt', 4, 1));
            S = cell2struct(S.ncst', 'points', 1);
            disp('Read shapefile with m_shaperead');
        end
        SG = [SG; S];
    end
end

%% Convert NaN-Delimited Vector to Struct
if isempty(SG)
    polygon_struct.outer = boubox;
    polygon_struct.inner = [];
    polygon_struct.mainland = [];
    polygon_struct.linestrings = [];
    return;
end

polygon_struct = struct('outer', boubox, 'inner', [], 'mainland', [], ...
                        'innerb', [], 'mainlandb', [], ...
                        'innerb_type', [], 'mainlandb_type', [], ...
                        'linestrings', []); % Store line strings
edges = Get_poly_edges(polygon_struct.outer);

%% Convert to Matrices and Separate Polygons from Line Strings
tmpC = struct2cell(SG)';
tmpC = tmpC(~cellfun(@isempty, tmpC(:,1)),:); % Remove empty polygons
tmpC = cellfun(@(row) row(:,1:2), tmpC, 'UniformOutput', false);

valid_polygons = {}; % Store valid polygons
line_strings = {};  % Store line strings

for i = 1:size(tmpC,1)
    points = tmpC{i,1};
    if size(points,2) > 2, points = points(:,1:2); end % Keep X, Y only

    % **Check for Polygon Closure Using Tolerance**
    first_point = points(1, :);
    distances = sqrt(sum((points(2:end, :) - first_point).^2, 2)); % Compute distances (excluding first)
    
    if any(distances < tolerance)  % **If any other point is within tolerance → It's a polygon**
        valid_polygons{end+1} = points; 
    else
        line_strings{end+1} = [points; NaN NaN]; % Store as a line string with NaN separator
    end
end

%% Store Line Strings in polygon_struct with NaN Separation
if ~isempty(line_strings)
    polygon_struct.linestrings = cell2mat(line_strings');
end

%% Classify Polygons into Mainland or Island
for points = valid_polygons
    points = points{1};
    area = shoelace(points(:,1), points(:,2)); % Compute polygon area
    inside_bbox = all(inpoly(points, polygon_struct.outer, edges));
    if inside_bbox && abs(area) >= 4 * h0^2
        polygon_struct.inner = [polygon_struct.inner; points; NaN NaN]; % Island
    elseif abs(area) >= 100 * h0^2
        polygon_struct.mainland = [polygon_struct.mainland; points; NaN NaN]; % Mainland
    end
end

%% Merge Overlapping Mainland & Inner Boundaries While Preserving NaNs
if exist('polyshape', 'file')
    if ~isempty(polygon_struct.mainland) && ~isempty(polygon_struct.inner)
        
        % Extract mainland polygons while preserving NaNs
        idx_m = find(isnan(polygon_struct.mainland(:,1)));
        idx_m = [0; idx_m; size(polygon_struct.mainland,1)+1]; 
        mainland_parts = arrayfun(@(j) polygon_struct.mainland(idx_m(j)+1:idx_m(j+1)-1,:), ...
                                  1:length(idx_m)-1, 'UniformOutput', false);
        
        % Extract inner polygons while preserving NaNs
        idx_i = find(isnan(polygon_struct.inner(:,1)));
        idx_i = [0; idx_i; size(polygon_struct.inner,1)+1]; 
        inner_parts = arrayfun(@(j) polygon_struct.inner(idx_i(j)+1:idx_i(j+1)-1,:), ...
                               1:length(idx_i)-1, 'UniformOutput', false);
        
        % Convert each part into polyshape and compute union
        warning('off', 'MATLAB:polyshape:repairedBySimplify');
        merged_polyshapes = [];
        for p1 = mainland_parts
            for p2 = inner_parts
                if ~isempty(p1{1}) && ~isempty(p2{1}) % Ensure non-empty parts
                    poly_m = polyshape(p1{1}(:,1), p1{1}(:,2));
                    poly_i = polyshape(p2{1}(:,1), p2{1}(:,2));
                    poly_merged = union(poly_m, poly_i);
                    merged_polyshapes = [merged_polyshapes; poly_merged.Vertices]; % Store merged polygons
                    merged_polyshapes = [merged_polyshapes; NaN NaN]; % Preserve NaN separators
                end
            end
        end
        % Restore warnings after execution
        warning('on', 'MATLAB:polyshape:repairedBySimplify');
        
        % Store back into polygon_struct.mainland
        polygon_struct.mainland = merged_polyshapes;
    end
end

%% Plot Results (Optional)
if plot_on >= 1
    figure(1); hold on;
    plot(polygon_struct.outer(:,1), polygon_struct.outer(:,2), 'k');
    plot(polygon_struct.inner(:,1), polygon_struct.inner(:,2), 'b');
    plot(polygon_struct.mainland(:,1), polygon_struct.mainland(:,2), 'r');
    if ~isempty(polygon_struct.linestrings)
        plot(polygon_struct.linestrings(:,1), polygon_struct.linestrings(:,2), 'c--'); % Line strings in cyan dashed
    end
    axis equal
end

end
