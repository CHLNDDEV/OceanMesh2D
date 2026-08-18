function polygon_struct = Read_shapefile(finputname, ~, bbox, h0, boubox, plot_on, ~, densify_outer)
%% ============================
%  Initialize Data Structures
% ============================
if nargin < 8 || isempty(densify_outer)
    densify_outer = 0;
end
SG = [];
loop = 1; minus = 0;
tolerance = 1e-5;
min_area_inner = 4 * h0^2;
min_area_mainland = 100 * h0^2;

if densify_outer
    % Densify the outer (bounding-box) polygon so gaps larger than half
    % the minimum edge length are filled in. Off by default: unnecessary
    % for most meshes and slower on large bounding boxes.
    [latout, lonout] = my_interpm(boubox(:,2), boubox(:,1), h0/2);
    boubox = [lonout, latout];
end

if bbox(1,2) > 180 && bbox(1,1) < 180, loop = 2; end
if all(bbox(1,:) > 180), minus = 1; end

%% ============================
%  Read Shapefile Data
% ============================
for fname = finputname
    for nn = 1:loop
        bboxt = bbox';
        if loop == 2, bboxt(2-((nn-1)*1),1) = 180 - (nn-1)*360; end
        if minus, bboxt(:,1) = bboxt(:,1) - 360; end

        try
            S = shaperead(fname{1}, 'BoundingBox', bboxt);
            S = rmfield(S, setdiff(fieldnames(S), {'X', 'Y'}));
        catch
            S = m_shaperead(fname{1}, reshape(bboxt', 4, 1));
            S = cell2struct(S.ncst', 'points', 1);
        end
        SG = [SG; S];
    end
end

%% ============================
%  Handle Empty Data
% ============================
if isempty(SG)
    polygon_struct = struct('outer', boubox, 'inner', [], 'mainland', [], ...
                            'linestrings', [], 'mainlandb', [], 'innerb', [], ...
                            'mainlandb_type', [], 'innerb_type', []);
    return;
end

polygon_struct = struct('outer', boubox, 'inner', [], 'mainland', [], ...
                        'linestrings', [], 'mainlandb', [], 'innerb', [], ...
                        'mainlandb_type', [], 'innerb_type', []);
edges = Get_poly_edges(polygon_struct.outer);

%% ============================
%  Extract & Classify Shapes
% ============================
tmpC = struct2cell(SG)';
tmpC = tmpC(~cellfun(@isempty, tmpC(:,1)),:);  
tmpC = cellfun(@(row) row(:,1:2), tmpC, 'UniformOutput', false);  

valid_polygons = {};
is_open_poly = false(0,1);
line_strings = {};

for i = 1:size(tmpC,1)
    points = tmpC{i,1};

    % **Check for Polygon Closure Using Tolerance**
    first_point = points(1, :);
    distances = sqrt(sum((points(2:end, :) - first_point).^2, 2));
    is_closed = any(distances < tolerance);

    if is_closed
        valid_polygons{end+1} = points; %#ok<AGROW>
        is_open_poly(end+1,1) = false; %#ok<AGROW>
    else
        line_strings{end+1} = [points; NaN NaN];  %#ok<AGROW>

        % Also carry not-exactly-closed features (e.g. a barrier
        % island/spit digitized as an open shoreline trace) through the
        % mainland/inner classification below as a pseudo-closed
        % polygon, so they aren't silently dropped from the mesh
        % boundary when linestrings aren't consumed (i.e. high_fidelity
        % isn't enabled). Mirrors the previous Read_shapefile's
        % behavior, which gave non-closed features a placeholder area
        % that always passed the size thresholds.
        valid_polygons{end+1} = points; %#ok<AGROW>
        is_open_poly(end+1,1) = true; %#ok<AGROW>
    end
end

if ~isempty(line_strings)
    polygon_struct.linestrings = cell2mat(line_strings');
end

%% ============================
%  Compute Polygon Areas and Assign to Inner or Mainland
% ============================
inner_polys = {}; 
mainland_polys = {}; 
innerb_types = {}; 
mainlandb_types = {}; 

% **Progress Bar for Large Datasets**
total_polygons = numel(valid_polygons);
if total_polygons > 1000
    f = waitbar(0, 'Processing polygons...');
end

for i = 1:total_polygons
    points = valid_polygons{i};
    if is_open_poly(i)
        % Not exactly closed: use a placeholder area, large enough to
        % clear min_area_inner/min_area_mainland regardless of h0, so
        % the feature is still classified as boundary geometry (see note
        % above where it was added to valid_polygons).
        area = 999;
    else
        area = shoelace(points(:,1), points(:,2));
    end
    inside_bbox = all(inpoly(points, polygon_struct.outer, edges));

    if inside_bbox && abs(area) >= min_area_inner
        inner_polys{end+1} = points;
        innerb_types{end+1} = 'inner';
    elseif abs(area) >= min_area_mainland
        mainland_polys{end+1} = points;
        mainlandb_types{end+1} = 'mainland';
    end

    % **Update Progress Bar**
    if exist('f', 'var') && mod(i, 500) == 0
        waitbar(i / total_polygons, f);
    end
end

if exist('f', 'var'), close(f); end  % Close progress bar if it exists

polygon_struct.inner = cell2mat(cellfun(@(x) [x; NaN NaN], inner_polys, 'UniformOutput', false)');
polygon_struct.mainland = cell2mat(cellfun(@(x) [x; NaN NaN], mainland_polys, 'UniformOutput', false)');

% Remove parts of inner and mainland overlapping with outer
polygon_struct.outer = [polygon_struct.outer; polygon_struct.mainland];

%% ============================
%  Initialize `mainlandb` and `innerb` Correctly
% ============================
polygon_struct.mainlandb = polygon_struct.mainland;  
polygon_struct.innerb = polygon_struct.inner;  
polygon_struct.mainlandb_type = mainlandb_types;
polygon_struct.innerb_type = innerb_types;

%% ============================
%  Plot Results (Optional)
% ============================
if plot_on >= 1
    figure(1); hold on;
    plot(polygon_struct.outer(:,1), polygon_struct.outer(:,2), 'k');
    plot(polygon_struct.inner(:,1), polygon_struct.inner(:,2), 'b');
    plot(polygon_struct.mainland(:,1), polygon_struct.mainland(:,2), 'r');
    if ~isempty(polygon_struct.linestrings)
        plot(polygon_struct.linestrings(:,1), polygon_struct.linestrings(:,2), 'c--');
    end
    axis equal
end

end
