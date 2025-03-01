function polygon_struct = Read_shapefile(finputname, ~, bbox, h0, boubox, plot_on, ~)

%% **User-Modifiable Parameters**
SMALL_POLYGON_FACTOR = 4;  % Factor for removing small polygons (e.g., 4 * h0^2)
LARGE_POLYGON_FACTOR = 100; % Factor for distinguishing mainland from small features
CLOSURE_TOLERANCE = 1e-2; % Tolerance for determining if a shape is a closed polygon

%% **Initialize Data Structures**
SG = []; % Store read polygons
loop = 1; minus = 0;
if bbox(1,2) > 180 && bbox(1,1) < 180, loop = 2; end
if all(bbox(1,:) > 180), minus = 1; end

%% **Read Shapefile Data**
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

%% **Handle Empty Case**
if isempty(SG)
    polygon_struct.outer = boubox;
    polygon_struct.inner = [];
    polygon_struct.mainland = [];
    polygon_struct.linestrings = [];
    return;
end

%% **Initialize polygon_struct**
polygon_struct = struct('outer', boubox, 'inner', [], 'mainland', [], ...
                        'innerb', [], 'mainlandb', [], ...
                        'innerb_type', [], 'mainlandb_type', [], ...
                        'linestrings', []); % Store line strings
edges = Get_poly_edges(polygon_struct.outer);

%% **Convert to Matrices and Separate Polygons from Line Strings**
tmpC = struct2cell(SG)';
tmpC = tmpC(~cellfun(@isempty, tmpC(:,1)),:); % Remove empty polygons

valid_polygons = {}; % Store valid polygons
line_strings = {};  % Store line strings

for i = 1:size(tmpC,1)
    points = tmpC{i,1};
    if size(points,2) > 2, points = points(:,1:2); end % Keep X, Y only

    % **Check if shape is a polygon (first and last point are within tolerance)**
    if norm(points(1,:) - points(end,:)) > CLOSURE_TOLERANCE
        line_strings{end+1} = [points; NaN NaN]; % Store as a line string with NaN separator
        continue; % Skip further processing
    end

    valid_polygons{end+1} = points; % Store as a valid polygon
end

%% **Store Line Strings in polygon_struct with NaN Separation**
if ~isempty(line_strings)
    polygon_struct.linestrings = cell2mat(line_strings');
end

%% **Classify Polygons into Mainland or Island**
for points = valid_polygons
    points = points{1};
    area = shoelace(points(:,1), points(:,2)); % Compute polygon area
    inside_bbox = all(inpoly(points, polygon_struct.outer, edges));

    % **Remove small polygons using SMALL_POLYGON_FACTOR**
    if inside_bbox && area >= SMALL_POLYGON_FACTOR * h0^2
        polygon_struct.inner = [polygon_struct.inner; points; NaN NaN]; % Island
    elseif area >= LARGE_POLYGON_FACTOR * h0^2
        polygon_struct.mainland = [polygon_struct.mainland; points; NaN NaN]; % Mainland
    end
end

%% **Merge Overlapping Mainland & Inner Boundaries**
if exist('polyshape', 'file')
    if ~isempty(polygon_struct.mainland) && ~isempty(polygon_struct.inner)
        poly_m = polyshape(polygon_struct.mainland);
        poly_i = polyshape(polygon_struct.inner);
        poly_merged = union(poly_m, poly_i);
        polygon_struct.mainland = poly_merged.Vertices;
    end
end

%% **Plot Results (Optional)**
if plot_on >= 1
    figure(1); hold on;
    plot(polygon_struct.outer(:,1), polygon_struct.outer(:,2), 'k');
    plot(polygon_struct.inner(:,1), polygon_struct.inner(:,2), 'b');
    plot(polygon_struct.mainland(:,1), polygon_struct.mainland(:,2), 'r');
    if ~isempty(polygon_struct.linestrings)
        plot(polygon_struct.linestrings(:,1), polygon_struct.linestrings(:,2), 'c--'); % Line strings in cyan dashed
    end
end

end
