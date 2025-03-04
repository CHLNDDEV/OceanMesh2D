classdef meshgen
    %   MESHGEN: Mesh generation class
    %   Handles input parameters to create a meshgen class object that can be
    %   used to build a msh class.
    %   Copyright (C) 2018  Keith Roberts & William Pringle
    %
    %   This program is free software: you can redistribute it and/or modify
    %   it under the terms of the GNU General Public License as published by
    %   the Free Software Foundation, either version 3 of the License, or
    %   (at your option) any later version.
    %
    %   This program is distributed in the hope that it will be useful,
    %   but WITHOUT ANY WARRANTY; without even the implied warranty of
    %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %   GNU General Public License for more details.
    %
    %   You should have received a copy of the GNU General Public License
    %   along with this program.  If not, see <http://www.gnu.org/licenses/>.

    properties
        % User-defined inputs and internal parameters
        fd            % Handle to distance function
        fh            % Handle to edge function(s)
        h0            % Minimum edge length (scalar or array)
        bbox          % Bounding box [xmin,ymin; xmax,ymax] (or cell array)
        pfix          % Fixed node positions (nfix x 2)
        egfix         % Fixed edge constraints (indices into pfix)
        fixboxes      % Flag indicating which boxes use fixed constraints
        plot_on       % Plotting flag (1 = on, 0 = off)
        nscreen       % Frequency of plotting and writing temp files
        bou           % Geodata object(s)
        ef            % Edgefx class instance(s)
        itmax         % Maximum number of iterations for mesh improvement
        outer         % Outer boundary (cell array if multiple boxes)
        inner         % Island boundaries (cell array)
        mainland      % Shoreline boundary (cell array)
        boubox        % Bounding box as a polygon (cell array)
        inpoly_flip   % Flag to flip the inpoly test for signed distances
        memory_gb     % Memory in GB for initial rejection method
        cleanup       % Flag/string to trigger mesh cleanup (default 'default')
        direc_smooth  % Flag to trigger direct smoothing during cleanup
        dj_cutoff     % Cutoff area fraction for deleting disjoint portions
        grd = msh();  % Mesh container (of type msh) for nodes p and triangles t
        qual          % Quality metrics for mesh elements
        qual_tol      % Tolerance for negligible mesh quality change
        proj          % Projection structure for m_map
        anno          % Approximate Nearest Neighbor search object(s)
        annData       % Data contained in the KD-tree(s)
        Fb            % Bathymetry data interpolant(s)
        enforceWeirs  % Flag to enforce weirs in mesh generation
        enforceMin    % Flag to enforce minimum edge length for all ef's
        improve_boundary  % Flag to improve the boundary representation
        high_fidelity     % Cell array or scalar flag for high-fidelity mesh
        delaunay_elim_on_exit % Flag to run delaunay_elim on exit
        improve_with_reduced_quality % Allow improvements with reduced quality
    end

    methods
        %% Constructor
        function obj = meshgen(varargin)
            % Parse input arguments and set defaults.
            p = inputParser;
            defval = 0; % placeholder for unspecified args
            addOptional(p,'h0',defval);
            addOptional(p,'bbox',defval);
            addOptional(p,'fh',defval);
            addOptional(p,'pfix',defval);
            addOptional(p,'egfix',defval);
            addOptional(p,'fixboxes',defval);
            addOptional(p,'inner',defval);
            addOptional(p,'outer',defval);
            addOptional(p,'mainland',defval);
            addOptional(p,'bou',defval);
            addOptional(p,'ef',defval);
            addOptional(p,'plot_on',defval);
            addOptional(p,'nscreen',defval);
            addOptional(p,'itmax',defval);
            addOptional(p,'memory_gb',1);
            addOptional(p,'cleanup',1);
            addOptional(p,'direc_smooth',1);
            addOptional(p,'dj_cutoff',0.25);
            addOptional(p,'big_mesh',defval);
            addOptional(p,'proj',defval);
            addOptional(p,'qual_tol',defval);
            addOptional(p,'enforceWeirs',0);
            addOptional(p,'enforceMin',1);
            addOptional(p,'delaunay_elim_on_exit',1);
            addOptional(p,'improve_with_reduced_quality',0);
            parse(p,varargin{:});
            inp = orderfields(p.Results,{'h0','bbox','enforceWeirs','enforceMin',...
                'delaunay_elim_on_exit','improve_with_reduced_quality',...
                'fh','inner','outer','mainland',...
                'bou','ef','egfix','pfix','fixboxes',...
                'plot_on','nscreen','itmax',...
                'memory_gb','qual_tol','cleanup',...
                'direc_smooth','dj_cutoff',...
                'big_mesh','proj'});

            % Loop through options and assign to object properties.
            fields = fieldnames(inp);
            for i = 1:numel(fields)
                switch fields{i}
                    case 'h0'
                        obj.h0 = inp.h0;
                    case 'fh'
                        if isa(inp.fh, 'function_handle')
                            obj.fh = inp.fh;
                        end
                    case 'bbox'
                        obj.bbox = inp.bbox;
                        if iscell(obj.bbox)
                            ob_min = obj.bbox{1}(:,1);
                            ob_max = obj.bbox{1}(:,2);
                            for ii = 2:length(obj.bbox)
                                if any(obj.bbox{ii}(:,1) < ob_min) || any(obj.bbox{ii}(:,2) > ob_max)
                                    error(['Outer bbox must contain all inner bboxes: inner box #' num2str(ii) ' violates this'])
                                end
                            end
                        end
                        if obj.bbox(1)==0, obj.bbox = []; end
                    case 'pfix'
                        obj.pfix = inp.pfix;
                        if ~isempty(obj.pfix) && obj.pfix(1) == 0, obj.pfix = []; end
                        if obj.enforceWeirs
                            for j = 1:length(obj.bou)
                                if ~isempty(obj.bou{j}.weirPfix)
                                    obj.pfix = [obj.pfix; obj.bou{j}.weirPfix];
                                end
                            end
                        end
                    case 'egfix'
                        obj.egfix = inp.egfix;
                        if ~isempty(obj.egfix) && obj.egfix(1)==0, obj.egfix = []; end
                        if obj.enforceWeirs
                            for j = 1:length(obj.bou)
                                if ~isempty(obj.bou{j}.weirEgfix) && ~isempty(obj.egfix)
                                    obj.egfix = [obj.egfix; obj.bou{j}.weirEgfix + max(obj.egfix(:))];
                                else
                                    obj.egfix = obj.bou{j}.weirEgfix;
                                end
                            end
                        end
                        obj.egfix = renumberEdges(obj.egfix);
                    case 'fixboxes'
                        obj.fixboxes = inp.fixboxes;
                    case 'bou'
                        if obj.outer ~= 0, continue; end
                        obj.outer = {}; obj.inner = {}; obj.mainland = {};
                        obj.bou = inp.bou;
                        if ~iscell(obj.bou)
                            obj.bou = {obj.bou};
                        end
                        for ee = 1:length(obj.bou)
                            arg = obj.bou{ee};
                            if isa(arg, 'geodata')
                                obj.high_fidelity{ee} = obj.bou{ee}.high_fidelity;
                                obj.outer{ee} = obj.bou{ee}.outer;
                                obj.inner{ee} = obj.bou{ee}.inner;
                                if ~isempty(obj.bou{ee}.Fb)
                                    obj.Fb{ee} = obj.bou{ee}.Fb;
                                end
                                if ~isempty(obj.inner{ee}) && obj.inner{ee}(1) ~= 0
                                    obj.outer{ee} = [obj.outer{ee}; obj.inner{ee}];
                                end
                                obj.mainland{ee} = obj.bou{ee}.mainland;
                                obj.boubox{ee} = obj.bou{ee}.boubox;
                                obj.inpoly_flip{ee} = obj.bou{ee}.inpoly_flip;
                            end
                        end
                    case 'ef'
                        tmp = inp.ef;
                        if isa(tmp, 'function_handle')
                            error('Please specify your edge function handle via the name/value pair fh');
                        end
                        obj.ef = tmp;
                        if ~iscell(obj.ef)
                            obj.ef = {obj.ef};
                        end
                        for ee = 1:length(obj.ef)
                            if isa(obj.ef{ee}, 'edgefx')
                                obj.bbox{ee} = obj.ef{ee}.bbox;
                            end
                        end
                        if iscell(obj.bbox)
                            ob_min = obj.bbox{1}(:,1);
                            ob_max = obj.bbox{1}(:,2);
                            for ii = 2:length(obj.bbox)
                                if any(obj.bbox{ii}(:,1) < ob_min) || any(obj.bbox{ii}(:,2) > ob_max)
                                    error(['Outer bbox must contain all inner bboxes: inner box #' num2str(ii) ' violates this'])
                                end
                            end
                        end
                        for ee = 1:length(obj.ef)
                            if isa(obj.ef{ee}, 'edgefx')
                                obj.h0(ee) = obj.ef{ee}.h0;
                            end
                        end
                        if length(obj.ef) > 1 && obj.enforceMin
                            obj.ef = enforce_min_ef(obj.ef);
                        end
                        obj.ef = smooth_outer(obj.ef, obj.Fb);
                        for ee = 1:length(obj.ef)
                            if isa(obj.ef{ee}, 'edgefx')
                                obj.fh{ee} = @(p) obj.ef{ee}.F(p);
                            end
                        end
                    case 'plot_on'
                        obj.plot_on = inp.plot_on;
                    case 'nscreen'
                        obj.nscreen = inp.nscreen;
                        if obj.nscreen ~= 0
                            obj.plot_on = 1;
                        else
                            obj.nscreen = 5;
                        end
                    case 'itmax'
                        obj.itmax = inp.itmax;
                        if obj.itmax == 0
                            obj.itmax = 100;
                            warning('No itmax specified; defaulting to 100');
                        end
                    case 'qual_tol'
                        obj.qual_tol = inp.qual_tol;
                        if obj.qual_tol == 0, obj.qual_tol = 0.01; end
                    case 'inner'
                        if ~isa(obj.bou, 'geodata')
                            obj.inner = inp.inner;
                        end
                    case 'outer'
                        if ~isa(obj.bou, 'geodata')
                            obj.outer = inp.outer;
                            if obj.inner(1) ~= 0
                                obj.outer = [obj.outer; obj.inner];
                            end
                        end
                    case 'mainland'
                        if ~isa(obj.bou, 'geodata')
                            obj.mainland = inp.mainland;
                        end
                    case 'memory_gb'
                        obj.memory_gb = inp.memory_gb;
                    case 'cleanup'
                        obj.cleanup = inp.cleanup;
                        if isempty(obj.cleanup) || obj.cleanup == 0
                            obj.cleanup = 'none';
                        elseif obj.cleanup == 1
                            obj.cleanup = 'default';
                        end
                    case 'dj_cutoff'
                        obj.dj_cutoff = inp.dj_cutoff;
                    case 'direc_smooth'
                        obj.direc_smooth = inp.direc_smooth;
                    case 'proj'
                        obj.proj = inp.proj;
                        if obj.proj == 0, obj.proj = 'equi'; end
                        if ~isempty(obj.bbox)
                            lon_mi = obj.bbox{1}(1,1) - obj.h0(1)/1110;
                            lon_ma = obj.bbox{1}(1,2) + obj.h0(1)/1110;
                            lat_mi = obj.bbox{1}(2,1) - obj.h0(1)/1110;
                            lat_ma = obj.bbox{1}(2,2) + obj.h0(1)/1110;
                        else
                            lon_mi = -180; lon_ma = 180; lat_mi = -90; lat_ma = 90;
                        end
                        dmy = msh();
                        dmy.p(:,1) = [lon_mi; lon_ma];
                        dmy.p(:,2) = [lat_mi; lat_ma];
                        setProj(dmy,1,obj.proj);
                    case 'enforceWeirs'
                        obj.enforceWeirs = inp.enforceWeirs;
                    case 'enforceMin'
                        obj.enforceMin = inp.enforceMin;
                    case 'delaunay_elim_on_exit'
                        obj.delaunay_elim_on_exit = inp.delaunay_elim_on_exit;
                    case 'improve_with_reduced_quality'
                        obj.improve_with_reduced_quality = inp.improve_with_reduced_quality;
                end
            end

            % Essential input checks.
            if any(obj.h0 == 0)
                error('h0 (minimum edge length) was not correctly specified!');
            end
            if isempty(obj.outer)
                error('No outer boundary specified!');
            end
            if isempty(obj.bbox)
                error('No bounding box specified!');
            end

            % Set default distance function and build ANN.
            obj.fd = @dpoly;
            obj = createANN(obj);
            global MAP_PROJECTION MAP_COORDS MAP_VAR_LIST
            obj.grd.proj    = MAP_PROJECTION;
            obj.grd.coord   = MAP_COORDS;
            obj.grd.mapvar  = MAP_VAR_LIST;

            % Generate breakline constraints from shoreline and island boundaries.
            [tpfix, tegfix] = obj.generateBreaklineConstraints();

            % Check if any entry in high_fidelity is 3
            if any(cellfun(@(x) isscalar(x) && x == 3, obj.high_fidelity))
                disp('     Pruning breakline connectivity due to high-fidelity mode 3...');

                % Set tolerance and angle threshold
                tol = (obj.h0(end)) / 111e3 ;
                angle_thresh = 30;

                % Apply CleanPSLG processing with timing
                disp('     Starting PSLG cleaning process...');

                % **Step 1: Merge Close Vertices**
                tic;
                pslg = CleanPSLG(tpfix, tegfix, tol, angle_thresh);
                pslg = pslg.mergeVertices();
                mergeTime = toc;
                fprintf('    Merged "too close" vertices in %.2f seconds\n', mergeTime);

                %figure;
                %pslg.plotPSLG();
                %title('After merge...');

                % **Step 2: Drop Intersecting Edges**
                tic;
                pslg = pslg.dropIntersectingEdges();
                dropTime = toc;
                fprintf('    Dropped intersecting edges in %.2f seconds\n', dropTime);
                % 
                % figure;
                % pslg.plotPSLG();
                % title('After drop intersect');

                % **Step 3: Prune Encroaching Edges**
                tic;
                pslg = pslg.pruneEncroachingEdges();
                pruneTime = toc;
                fprintf('    Pruned encroaching edges in %.2f seconds\n', pruneTime);

                % figure;
                % pslg.plotPSLG();
                % title('After prune encroach...');

                % **Total time**
                totalTime = mergeTime + dropTime + pruneTime;
                fprintf('   Total PSLG processing time: %.2f seconds\n', totalTime);

                % Update points and edges
                tpfix = pslg.Vertices;
                tegfix = pslg.Segments;
            end

            % Check for duplicate fixed points.
            if ~isempty(tpfix)
                checkfixp = setdiff(tpfix, fixmesh(tpfix), 'rows');
                if ~isempty(checkfixp)
                    error('Duplicate fixed points detected, cannot proceed');
                end
            end

            % Update object with user-supplied and generated constraints.
            obj.pfix = [obj.pfix; tpfix];
            if isempty(obj.egfix)
                obj.egfix = tegfix;
            else
                obj.egfix = [obj.egfix; tegfix + max(obj.egfix(:))];
            end

            if ~isempty(obj.pfix)
                disp(['Using ', num2str(size(obj.pfix,1)), ' fixed points.']);
            end
            if ~isempty(obj.egfix)
                if max(obj.egfix(:)) > size(obj.pfix,1)
                    error('FATAL: egfix indices exceed the number of fixed points.');
                end
                disp(['Using ', num2str(size(obj.egfix,1)), ' fixed edges.']);
                warning('Please verify fixed constraints using plot.');
            end
        end



        %% Private method: generateBreaklineConstraints
        function [tpfix, tegfix] = generateBreaklineConstraints(obj)
            % Generates breakline constraints (fixed points and edges) from the
            % provided mainland, inner (island) boundaries, and line strings.
            %
            % Outputs:
            %   tpfix  - Accumulated fixed points from polygon boundaries and line strings.
            %   tegfix - Corresponding fixed edge constraints.

            tpfix = [];
            tegfix = [];

            for box_num = 1:length(obj.h0)
                % Collect polygon data from mainland and inner boundaries.
                polys = {};
                if ~isempty(obj.mainland{box_num})
                    polys{end+1} = obj.mainland{box_num};
                end
                if ~isempty(obj.inner{box_num}) && obj.inner{box_num}(1) ~= 0
                    polys{end+1} = obj.inner{box_num};
                end

                % Collect NaN-separated line strings
                lineStrings = {};
                if ~isempty(obj.bou{box_num}.linestrings)
                    lineStrings{end+1} = obj.bou{box_num}.linestrings; % Store line strings
                end

                % Process only if high-fidelity is enabled for this box.
                if obj.high_fidelity{box_num}
                    if isequal(obj.cleanup, 1)
                        warning('Disabling cleanup since high-fidelity mode is active.');
                        obj.cleanup = 0;
                    end
                    disp(['Redistributing vertices for box #', num2str(box_num)]);

                    % Process Polygons: Concatenate and remove NaNs
                    poly_all = cell2mat(polys');
                    D = nandelim_to_cell(poly_all);  % Split polygons at NaNs.
                    D = D(~cellfun(@(p) all(isnan(p(:))), D));  % Remove empty segments
                    areas = cellfun(@(d) polyarea(d(:,1), d(:,2)), D);
                    [~, uniqueIdx] = uniquetol(areas, 1e-12);
                    D = D(uniqueIdx);

                    % Process each polygon separately
                    for i = 1:length(D)
                        current_poly = D{i};
                        if size(current_poly,1) > 2
                            polyEdges = Get_line_edges(current_poly);
                            [pts, bnde] = filter_polygon_constraints(current_poly, polyEdges, obj.boubox, box_num);
                            if isempty(bnde), continue; end

                            % Split polygon if necessary.
                            polySplits = extdom_polygon(bnde, pts, 0, 1);
                            for j = 1:length(polySplits)
                                points = unique(polySplits{j}, 'rows', 'stable');

                                % Generate fixed constraints based on high_fidelity option.
                                if obj.high_fidelity{box_num} == 2
                                    tmp_pfix = points;
                                    tmp_egfix = Get_poly_edges([points; NaN, NaN]);
                                elseif obj.high_fidelity{box_num} == 1 || obj.high_fidelity{box_num} == 3
                                    [tmp_pfix, tmp_egfix] = mesh1d(points, obj.fh, obj.h0./111e3, [], obj.boubox, box_num, []);
                                else
                                    continue;
                                end

                                if size(tmp_pfix,1) > 2
                                    [tmp_pfix, tmp_egfix] = fixgeo2(tmp_pfix, tmp_egfix);
                                    if max(tmp_egfix(:)) ~= size(tmp_pfix,1), continue; end
                                    tpfix = [tpfix; tmp_pfix];
                                    if isempty(tegfix)
                                        tegfix = tmp_egfix;
                                    else
                                        tegfix = [tegfix; tmp_egfix + max(tegfix(:))];
                                    end
                                end
                            end
                        end
                    end

                    % Process NaN-separated Line Strings
                    for i = 1:length(lineStrings)
                        current_lines = lineStrings{i};
                        if size(current_lines,1) > 1
                            % Convert NaN-separated line strings into individual segments
                            lineSegments = nandelim_to_cell(current_lines); % Split at NaNs
                            lineSegments = lineSegments(~cellfun(@isempty, lineSegments)); % Remove empty segments

                            for j = 1:length(lineSegments)
                                lineSeg = unique(lineSegments{j}, 'rows', 'stable');
                                if size(lineSeg,1) < 2, continue; end % Ensure it's valid

                                % Generate fixed constraints
                                if obj.high_fidelity{box_num} == 2
                                    tmp_pfix = lineSeg;
                                    tmp_egfix = Get_poly_edges([lineSeg; NaN, NaN]);
                                elseif obj.high_fidelity{box_num} == 1 || obj.high_fidelity{box_num} == 3
                                    [tmp_pfix, tmp_egfix] = mesh1d(lineSeg, obj.fh, (3*obj.h0)./111e3, [], obj.boubox, box_num, []);
                                else
                                    continue;
                                end

                                if size(tmp_pfix,1) > 1
                                    [tmp_pfix, tmp_egfix] = fixgeo2(tmp_pfix, tmp_egfix);
                                    if max(tmp_egfix(:)) ~= size(tmp_pfix,1), continue; end
                                    tpfix = [tpfix; tmp_pfix];
                                    if isempty(tegfix)
                                        tegfix = tmp_egfix;
                                    else
                                        tegfix = [tegfix; tmp_egfix + max(tegfix(:))];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % Final adjustment: re-index and remove duplicate fixed constraints.
            [tpfix, tegfix] = fixgeo2(tpfix, tegfix);
        end


        function obj = plot(obj)
            if ~isempty(obj.pfix) && ~isempty(obj.egfix)
                figure; hold on; grid on;

                % Set colors and styles
                bboxColor = [0, 0.7, 0];    % Green for bounding boxes
                outerColor = [0, 0, 1];      % Blue for outer boundaries
                innerColor = [1, 0, 0];      % Red for inner boundaries
                edgeColor = [0, 0, 0];       % Black for fixed edges
                fixedPointColor = [0, 0, 0]; % Black for fixed points

                % Loop over bounding boxes and plot constraints
                for box_number = 1:length(obj.boubox)
                    iboubox = obj.boubox{box_number};

                    % Plot bounding box with slight transparency
                    plot(iboubox(:,1), iboubox(:,2), '-', 'Color', bboxColor, 'LineWidth', 2, ...
                        'DisplayName', 'Bounding Box', 'HandleVisibility', 'off');

                    % Plot outer boundary (if available)
                    touter = obj.outer(box_number);
                    if ~isempty(touter)
                        tedges = Get_poly_edges(touter{1});
                        [touter, ~] = filter_polygon_constraints(touter{1}, tedges, obj.boubox, box_number);
                        plot(touter(:,1), touter(:,2), 'o', 'Color', outerColor, 'MarkerSize', 5, ...
                            'DisplayName', 'Outer Boundary');
                    end

                    % Plot inner boundary (if available)
                    tinner = obj.inner(box_number);
                    if ~isempty(tinner) && ~isempty(tinner{1})
                        tedges = Get_poly_edges(tinner{1});
                        [tinner, ~] = filter_polygon_constraints(tinner{1}, tedges, obj.boubox, box_number);
                        plot(tinner(:,1), tinner(:,2), 'x', 'Color', innerColor, 'MarkerSize', 5, ...
                            'DisplayName', 'Inner Boundary');
                    end
                end

                % Plot fixed edges and points
                if exist('drawedge2', 'file') == 2
                    drawedge2(obj.pfix, obj.egfix, edgeColor);
                else
                    % Plot fixed points
                    scatter(obj.pfix(:,1), obj.pfix(:,2), 40, fixedPointColor, 'filled', ...
                        'DisplayName', 'Fixed Points');

                    % Plot fixed edges with a thicker line
                    for i = 1:size(obj.egfix, 1)
                        plot(obj.pfix(obj.egfix(i,:), 1), obj.pfix(obj.egfix(i,:), 2), '-', ...
                            'Color', edgeColor, 'LineWidth', 2.5, 'DisplayName', 'Fixed Edges');
                    end
                end

                % Improve plot formatting
                axis equal;
                title('Constrained Breaklines and Mesh Constraints', 'FontWeight', 'bold');
                xlabel('Longitude', 'FontSize', 12);
                ylabel('Latitude', 'FontSize', 12);

                % Remove duplicate legend entries
                legendEntries = findobj(gca, '-property', 'DisplayName');
            else
                disp('No constraints to plot!');
            end
        end


        function obj = createANN(obj)
            box_vec = 1:length(obj.bbox);
            for box_num = box_vec
                if ~iscell(obj.outer)
                    dataset = obj.outer;
                    dataset(isnan(obj.outer(:,1)),:) = [];
                else
                    dataset = obj.outer{box_num};
                    dataset(isnan(obj.outer{box_num}(:,1)),:) = [];
                end
                if all(abs(obj.bbox{box_num}(1,:)) == 180)
                    dataset(abs(dataset(:,1)) > 180-1e-6,:) = [];
                    dataset(abs(dataset(:,1)) < 1e-6,:) = [];
                end
                [dataset(:,1),dataset(:,2)] = m_ll2xy(dataset(:,1),dataset(:,2));
                dataset(isnan(dataset(:,1)),:) = [];
                dmy = ann(dataset');
                obj.anno{box_num} = dmy;
                obj.annData{box_num} = dataset;
            end
        end

        function mesh_out = collapse_thin_triangles(obj, aspect_ratio_threshold)
            % Identify and collapse thin triangles in the mesh
            % aspect_ratio_threshold: Defines what is considered "thin"

            tri = obj.t; % Get triangle connectivity
            nodes = obj.p; % Get node coordinates

            num_tri = size(tri, 1);

            for i = 1:num_tri
                % Get triangle node indices
                n1 = tri(i, 1);
                n2 = tri(i, 2);
                n3 = tri(i, 3);

                % Compute edge lengths
                e1 = norm(nodes(n2, :) - nodes(n1, :)); % Edge 1-2
                e2 = norm(nodes(n3, :) - nodes(n2, :)); % Edge 2-3
                e3 = norm(nodes(n1, :) - nodes(n3, :)); % Edge 3-1

                % Find longest edge
                [longest_edge, idx] = max([e1, e2, e3]);

                % Compute aspect ratio (shortest to longest)
                shortest_edge = min([e1, e2, e3]);
                aspect_ratio = shortest_edge / longest_edge;

                % If the triangle is too thin, collapse it
                if aspect_ratio < aspect_ratio_threshold
                    % Find the node opposite the longest edge
                    switch idx
                        case 1, opposite_node = n3; edge_nodes = [n1, n2];
                        case 2, opposite_node = n1; edge_nodes = [n2, n3];
                        case 3, opposite_node = n2; edge_nodes = [n3, n1];
                    end

                    % Compute midpoint of longest edge
                    midpoint = mean(nodes(edge_nodes, :), 1);

                    % Move the opposite node to the midpoint (collapsing the triangle)
                    nodes(opposite_node, :) = midpoint;
                end
            end

            % Update mesh with modified nodes
            mesh_out = mesh;
            mesh_out.p = nodes;
        end


        function  obj = build(obj)
            % 2-D Mesh Generator using Distance Functions.
            % Checking existence of major inputs
            %%
            warning('off','all')
            %%
            tic
            it = 1 ;
            Re = 6378.137e3;
            geps = 1e-12*min(obj.h0)/Re;
            deps = sqrt(eps);
            ttol=0.1; Fscale = 1.2; deltat = 0.1;
            delIT = 0 ; delImp = 2;
            imp = 10; % number of iterations to do mesh improvements (delete/add)
            EXIT_QUALITY = 0.30; % minimum quality necessary to terminate if iter < itmax

            % unpack initial points.
            p = obj.grd.p;
            if isempty(p)
                disp('Forming initial point distribution...');
                % loop over number of boxes
                for box_num = 1:length(obj.h0)
                    disp(['    for box #' num2str(box_num)]);
                    % checking if cell or not and applying local values
                    h0_l = obj.h0(box_num);
                    max_r0 = 1/h0_l^2;
                    if ~iscell(obj.bbox)
                        bbox_l = obj.bbox'; % <--we must tranpose this!
                    else
                        bbox_l = obj.bbox{box_num}'; % <--tranpose!
                    end
                    if ~iscell(obj.fh)
                        fh_l = obj.fh;
                    else
                        fh_l = obj.fh{box_num};
                    end
                    % Lets estimate the num_points the distribution will be
                    num_points = ceil(2/sqrt(3)*prod(abs(diff(bbox_l)))...
                        /(h0_l/111e3)^2);
                    noblks = ceil(num_points*2*8/obj.memory_gb*1e-9);
                    len = abs(bbox_l(1,1)-bbox_l(2,1));
                    blklen = len/noblks;
                    st = bbox_l(1,1) ; ed = st + blklen; ns = 1;
                    %% 1. Create initial distribution in bounding box
                    %% (equilateral triangles)
                    for blk = 1 : noblks
                        if blk == noblks
                            ed = bbox_l(2,1);
                        end
                        ys = bbox_l(1,2);
                        ny = floor(1e3*m_lldist(repmat(0.5*(st+ed),2,1),...
                            [ys;bbox_l(2,2)])/h0_l);
                        dy = diff(bbox_l(:,2))/ny;
                        ns = 1;
                        % start at lower left and make grid going up to
                        % north latitude
                        for ii = 1:ny+1
                            if st*ed < 0
                                nx = floor(1e3*m_lldist([st;0],...
                                    [ys;ys])/(2/sqrt(3)*h0_l)) + ...
                                    floor(1e3*m_lldist([0;ed],...
                                    [ys;ys])/(2/sqrt(3)*h0_l));

                            else
                                nx = floor(1e3*m_lldist([st;ed],...
                                    [ys;ys])/(2/sqrt(3)*h0_l));
                            end
                            ne = ns+nx-1;
                            if mod(ii,2) == 0
                                % no offset
                                x(ns:ne) = linspace(st,ed,nx);
                            else
                                % offset
                                dx = (ed-st)/nx;
                                x(ns:ne) = linspace(st+0.5*dx,ed,nx);
                            end
                            % tolerance
                            if ii == (ny + 1)
                                y(ns:ne) = ys - eps;
                            else
                                y(ns:ne) = ys;
                            end
                            ns = ne+1; ys = ys + dy;

                        end

                        st = ed;
                        ed = st + blklen;
                        p1 = [x(:) y(:)]; clear x y


                        %% 2. Remove points outside the region, apply the rejection method
                        p1 = p1(feval(obj.fd,p1,obj,box_num) < geps,:);     % Keep only d<0 points
                        r0 = 1./feval(fh_l,p1).^2;                          % Probability to keep point
                        p1 = p1(rand(size(p1,1),1) < r0/max_r0,:);          % Rejection method
                        p  = [p; p1];                                       % Adding p1 to p
                    end
                    if box_num == 1
                        % add points along the outermost polygon to fill
                        % outer extent more quickly.
                        outer_temp = obj.outer{1};
                        Inan = find(isnan(outer_temp(:,1)),1,'first');
                        p1 = outer_temp(1:Inan-1,:);
                        p1 = p1(feval(obj.fd,p1,obj,1) < geps,:);     % Keep only d<0 points
                        r0 = 1./feval(fh_l, p1).^2;                         % Probability to keep point
                        p1 = p1(rand(size(p1,1),1) < r0/max_r0,:);          % Rejection method
                        p = [p; p1];                                        % Adding p1 to p
                    end
                end
            else
                disp('User-supplied initial points!');
                obj.grd.b = [];
                h0_l = obj.h0(end); % finest h0 (in case of a restart of meshgen.build).
            end

            nfix = length(obj.pfix); negfix = length(obj.egfix);
            if ~isempty(obj.pfix); p = [obj.pfix; p]; end
            % kjr July 2023, set to these values for better convg.
            if nfix > 0
                Fscale=1.1;
                deltat=0.10;
            end

            % Check if any boxes are set to high-fidelity
            % If so turn off heal_fixed_edges
            HIGH_FIDELITY_MODE = 0;
            for i = 1 : length(obj.h0)
                if obj.high_fidelity{i}
                    HIGH_FIDELITY_MODE = 1;
                end
            end

            N = size(p,1); % Number of points N
            disp(['Number of initial points after rejection is ',num2str(N)]);
            %% Iterate
            pold = inf;                                                    % For first iteration
            if obj.plot_on >= 1
                clf,view(2),axis equal;
            end
            toc
            fprintf(1,' ------------------------------------------------------->\n') ;
            disp('Begin iterating...');
            while 1
                tic
                if ~mod(it,obj.nscreen) && delIT == 0
                    disp(['Iteration = ' num2str(it)]) ;
                end
                % 3. Retriangulation by the Delaunay algorithm
                if max(sqrt(sum((p(1:size(pold,1),:)-pold).^2,2))/h0_l*111e3) > ttol         % Any large movement?
                    if it > 1
                        p = fixmesh([obj.pfix; p]);
                    else
                        p = fixmesh(p);                                        % Ensure only unique points.
                    end
                    N = size(p,1); pold = p;                               % Save current positions
                    [t,p] = delaunay_elim(p,obj.fd,geps,0);                % Delaunay with elimination
                    if isempty(t)
                        disp('Exiting')
                        return
                    end
                    % Getting element quality and check "goodness"
                    if exist('pt','var'); clear pt; end
                    [pt(:,1),pt(:,2)] = m_ll2xy(p(:,1),p(:,2));
                    tq = gettrimeshquan( pt, t);
                    mq_m = mean(tq.qm);
                    mq_l = min(tq.qm);
                    mq_s = std(tq.qm);
                    mq_l3sig = mq_m - 3*mq_s;
                    obj.qual(it,:) = [mq_m,mq_l3sig,mq_l];

                    % If not allowing improvements with reduction in quality
                    % ..or..
                    % If not allowing improvements with reduction in quality
                    % then if the number of points significantly decreased
                    % due to a mesh improvement iteration, then rewind.
                    if ~mod(it,imp+1) && ((obj.qual(it,1) - obj.qual(it-1,1) < -0.10)  || ...
                            (~obj.improve_with_reduced_quality && ...
                            (N - length(p_before_improve))/length(p_before_improve) < -0.10))

                        disp('Mesh improvement was unsuccessful...rewinding...');
                        p = p_before_improve;
                        N = size(p,1);                                     % Number of points changed
                        pold = inf;
                        it = it + 1;
                        continue
                    else
                        N = size(p,1); pold = p;                           % Assign number of points and save current positions
                    end
                    % 4. Describe each bar by a unique pair of nodes.
                    bars = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];           % Interior bars duplicated
                    bars = unique(sort(bars,2),'rows');                    % Bars as node pairs
                    % 5. Graphical output of the current mesh
                    if obj.plot_on >= 1 && (mod(it,obj.nscreen)==0 || it == 1)
                        cla,m_triplot(p(:,1),p(:,2),t)
                        m_grid
                        title(['Iteration = ',num2str(it)]);
                        if negfix > 0
                            m_plot(reshape(obj.pfix(obj.egfix,1),[],2)',...
                                reshape(obj.pfix(obj.egfix,2),[],2)','r-')
                        end
                        if nfix > 0
                            m_plot(obj.pfix(:,1),obj.pfix(:,2),'b.')
                        end
                        hold on ;
                        axis manual
                        drawnow
                    end
                end
                % Getting element quality and check goodness
                if exist('pt','var'); clear pt; end
                [pt(:,1),pt(:,2)] = m_ll2xy(p(:,1),p(:,2));
                tq = gettrimeshquan( pt, t);
                mq_m = mean(tq.qm);
                mq_l = min(tq.qm);
                mq_s = std(tq.qm);
                mq_l3sig = mq_m - 3*mq_s;
                obj.qual(it,:) = [mq_m,mq_l3sig,mq_l];

                % Improve the quality of triangles next to fixed edges by
                % deleting the point part of thin triangles without the fixed
                % point in it. Thin triangles have poor geometric quality <
                % 10%.
                if ~isempty(obj.egfix) && ~mod(it,delImp) && ~HIGH_FIDELITY_MODE
                    del = heal_fixed_edges(p,t,obj.egfix) ;
                    if ~isempty(del)
                        delIT = delIT + 1 ;
                        if delIT < 5
                            p(del,:)= [];
                            pold = inf;
                            disp(['Deleting ',num2str(length(del)),...
                                ' points close to fixed edges']);
                            continue;
                        else
                            % Abandon strategy..if it will not terminate
                            disp('Moving to next iteration');
                        end
                    end
                    delIT = 0 ;
                end

                % Termination quality, mesh quality reached is copacetic.
                qual_diff = mq_l3sig - obj.qual(max(1,it-imp),2);
                if ~mod(it,imp)
                    if mq_l > EXIT_QUALITY
                        % Do the final elimination of small connectivity
                        if obj.delaunay_elim_on_exit
                        end
                        disp('Quality of mesh is good enough, exit')
                        close all;
                        break;
                    end
                end
                % Saving a temp mesh
                if ~mod(it,obj.nscreen) && delIT == 0
                    disp(['Number of nodes is ' num2str(length(p))])
                    disp(['Mean mesh quality is ' num2str(mq_m)])
                    disp(['Min mesh quality is ' num2str(mq_l)])
                    disp(['3rd sigma lower mesh quality is ' num2str(mq_l3sig)])
                    tempp = p; tempt = t;
                    save('Temp_grid.mat','it','tempp','tempt');
                    clearvars tempp tempt
                end
                % 6. Move mesh points based on bar lengths L and forces F
                barvec = pt(bars(:,1),:)- pt(bars(:,2),:);                 % List of bar vectors
                if strcmp(obj.grd.proj.name,'UTM')
                    % UTM is already in meters (useful for small domains)
                    L = sqrt(sum(barvec.^2,2))*Re;
                else
                    % Get spherical earth distances
                    long   = zeros(length(bars)*2,1);
                    lat    = zeros(length(bars)*2,1);
                    long(1:2:end) = p(bars(:,1),1);
                    long(2:2:end) = p(bars(:,2),1);
                    lat(1:2:end)  = p(bars(:,1),2);
                    lat(2:2:end)  = p(bars(:,2),2);
                    L = m_lldist(long,lat); L = L(1:2:end)*1e3;            % L = Bar lengths in meters
                end
                ideal_bars = 0.5*(pt(bars(:,1),:) + pt(bars(:,2),:));      % Used to determine what bars are in bbox
                [ideal_bars(:,1),ideal_bars(:,2)] = ...                    % needs to be in non-projected
                    m_xy2ll(ideal_bars(:,1),ideal_bars(:,2));              % coordinates
                hbars = 0*ideal_bars(:,1);


                for box_num = 1:length(obj.h0)                             % For each bbox, find the bars that are in and calculate
                    if ~iscell(obj.fh)                                     % their ideal lengths.
                        fh_l = obj.fh;
                    else
                        fh_l = obj.fh{box_num};
                    end
                    h0_l = obj.h0(box_num);
                    if box_num > 1
                        h0_l = h0_l/111e3;                                 % create buffer to evalulate fh between nests
                        iboubox = obj.boubox{box_num}(1:end-1,:) ;
                        inside = inpoly(ideal_bars,iboubox) ;
                    else
                        inside = true(size(hbars));
                    end
                    hbars(inside) = feval(fh_l,ideal_bars(inside,:));      % Ideal lengths
                end


                L0 = hbars*Fscale*median(L)/median(hbars);                 % L0 = Desired lengths using ratio of medians scale factor
                LN = L./L0;                                                % LN = Normalized bar lengths


                % Mesh improvements (deleting and addition)
                p_before_improve = p;
                if ~mod(it,imp) %
                    nn = []; pst = [];
                    if abs(qual_diff) < imp*obj.qual_tol && ...
                            (obj.improve_with_reduced_quality || qual_diff > 0)

                        % Remove elements with small connectivity
                        nn = get_small_connectivity(p,t);
                        disp(['Deleting ' num2str(length(nn)) ' due to small connectivity'])


                        % Remove points that are too close (< LN = 0.5)
                        if any(LN < 0.5)
                            % Do not delete pfix too close.
                            nn1 = setdiff(reshape(bars(LN < 0.5,:),[],1),[(1:nfix)']);
                            disp(['Deleting ' num2str(length(nn1)) ' points too close together'])
                            nn = unique([nn; nn1]);
                        end


                        % Split long edges however many times to
                        % better lead to LN of 1
                        if any(LN > 2)
                            nsplit = floor(LN);
                            nsplit(nsplit < 1) = 1;
                            adding = 0;
                            % Split once
                            for jj = 2:2
                                il = find(nsplit >= jj);
                                xadd = zeros(length(il),jj-1);
                                yadd = zeros(length(il),jj-1);
                                for jjj = 1 : length(il)
                                    deltax = (p(bars(il(jjj),2),1)- p(bars(il(jjj),1),1))/jj;
                                    deltay = (p(bars(il(jjj),2),2)- p(bars(il(jjj),1),2))/jj;
                                    xadd(jjj,:) = p(bars(il(jjj),1),1) + (1:jj-1)*deltax;
                                    yadd(jjj,:) = p(bars(il(jjj),1),2) + (1:jj-1)*deltay;
                                end
                                pst = [pst; xadd(:) yadd(:)];
                                adding = numel(xadd) + adding;
                            end
                            disp(['Adding ',num2str(adding) ,' points.'])
                        end
                    end
                    if ~isempty(nn) || ~isempty(pst)
                        % Doing the actual subtracting and add
                        p(nn,:)= [];
                        p = [p; pst];
                        pold = inf;
                        it = it + 1;
                        continue;
                    end
                end


                F    = (1-LN.^4).*exp(-LN.^4)./LN;                         % Bessens-Heckbert edge force
                F(isinf(F)) = 0;
                Fvec = F*[1,1].*barvec;


                Ftot = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
                Ftot(1:nfix,:) = 0;                                        % Force = 0 at fixed points

                pt = pt + deltat*Ftot;                                     % Update node positions


                [p(:,1),p(:,2)] = m_xy2ll(pt(:,1),pt(:,2));


                %7. Bring outside points back to the boundary
                d = feval(obj.fd,p,obj,[],1); ix = d > 0;                  % Find points outside (d>0)
                ix(1:nfix) = 0;
                alpha = 1.0;
                for ib = 1 : 1 % length(obj.improve_boundary)
                    if sum(ix) > 0
                        pn = p(ix,:) + deps;
                        dgradx = (feval(obj.fd,[pn(:,1),p(ix,2)],obj,[])...%,1)...
                            -d(ix))/deps; % Numerical
                        dgrady = (feval(obj.fd,[p(ix,1),pn(:,2)],obj,[])...%,1)...
                            -d(ix))/deps; % gradient
                        dgrad2 = dgradx.^+2 + dgrady.^+2;
                        dgrad2(dgrad2 < eps) = eps;
                        p(ix,:) = p(ix,:) - alpha*[d(ix).*dgradx./dgrad2,...
                            d(ix).*dgrady./dgrad2];
                    end
                    alpha = alpha / 0.5;
                end


                % 8. Termination criterion: Exceed itmax
                it = it + 1 ;


                if ( it > obj.itmax )
                    % Do the final deletion of small connectivity
                    if obj.delaunay_elim_on_exit
                    end
                    disp('too many iterations, exit')
                    close all;
                    break ;
                end
                toc
            end
            %%
            warning('on','all')
            %%
            disp('Finished iterating...');
            fprintf(1,' ------------------------------------------------------->\n') ;


            %% Doing the final cleaning and fixing to the mesh...
            % Always save the mesh!
            save('Precleaned_grid.mat','it','p','t');

            % Clean up the mesh if specified
            if ~strcmp(obj.cleanup,'none')
                % Put the mesh class into the grd part of meshgen and clean
                obj.grd.p = p; obj.grd.t = t;
                [obj.grd,qout] = clean(obj.grd,obj.cleanup,...
                    'nscreen',obj.nscreen,'djc',obj.dj_cutoff,...
                    'pfix',obj.pfix);
                obj.grd.pfix = obj.pfix ;
                obj.grd.egfix= obj.egfix ;
                obj.grd.egfix= obj.egfix ;
                obj.qual(end+1,:) = qout;
            else
                % Fix mesh on the projected space
                [p(:,1),p(:,2)] = m_ll2xy(p(:,1),p(:,2));
                [p,t] = fixmesh(p,t);
                [p(:,1),p(:,2)] = m_xy2ll(p(:,1),p(:,2));
                % Put the mesh class into the grd part of meshgen
                obj.grd.p = p; obj.grd.t = t;
                obj.grd.pfix = obj.pfix ;
                obj.grd.egfix= obj.egfix ;
            end


            % Check element order, important for the global meshes crossing
            % -180/180 boundary
            obj.grd = CheckElementOrder(obj.grd);


            if obj.plot_on
                figure; plot(obj.qual,'linewi',2);
                hold on
                % plot the line dividing cleanup and distmesh
                plot([it it],[0 1],'--k')
                xticks(1:5:obj.itmax);
                xlabel('Iterations'); ylabel('Geometric element quality');
                title('Geometric element quality with iterations');
                set(gca,'FontSize',14);
                legend('q_{mean}','q_{mean}-q_{3\sigma}', 'q_{min}','Location','best');
                grid minor
            end
            return;
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Auxiliary subfunctions %
            %%%%%%%%%%%%%%%%%%%%%%%%%%


            function [t,p] = delaunay_elim(p,fd,geps,final)
                % Removing mean to reduce the magnitude of the points to
                % help the convex calc
                if exist('pt1','var'); clear pt1; end
                [pt1(:,1),pt1(:,2)] = m_ll2xy(p(:,1),p(:,2));
                if isempty(obj.egfix)
                    p_s  = pt1 - repmat(mean(pt1),[N,1]);
                    TR   = delaunayTriangulation(p_s);
                else
                    TR   = delaunayTriangulation(pt1(:,1),pt1(:,2),obj.egfix);
                    pt1  = TR.Points;
                end
                for kk = 1:final+1
                    if kk > 1
                        % Perform the following below upon exit from the mesh
                        % generation algorithm
                        nn = get_small_connectivity(pt1,t);
                        nn1 = [];
                        nn = unique([nn; nn1]) ;
                        TR.Points(nn,:) = [];
                        pt1(nn,:) = [];
                    end
                    t = TR.ConnectivityList;
                    pmid = squeeze(mean(reshape(pt1(t,:),[],3,2),2));      % Compute centroids
                    [pmid(:,1),pmid(:,2)] = m_xy2ll(pmid(:,1),pmid(:,2));  % Change back to lat lon
                    t    = t(feval(fd,pmid,obj,[]) < -geps,:);             % Keep interior trianglesi
                end
                if length(pt1) ~= length(p)
                    clear p
                    [p(:,1),p(:,2)] = m_xy2ll(pt1(:,1),pt1(:,2));
                end
            end


            function nn = get_small_connectivity(p,t)
                % Get node connectivity (look for 4)
                [~, enum] = VertToEle(t);
                % Make sure they are not boundary nodes
                bdbars = extdom_edges2(t, p);
                bdnodes = unique(bdbars(:));
                I = find(enum <= 4);
                nn = setdiff(I',[(1:nfix)';bdnodes]);                      % and don't destroy pfix or egfix!
                return;
            end

            function del = heal_fixed_edges(p,t,egfix)
                % kjr april2019
                % if there's a triangle with a low geometric quality that
                % contains a fixed edge, remove the non-fixed vertex
                % perform this on every other iteration to allow non-fixed
                % points to create equilateral triangles nearby the locked
                % edge.
                % returns points IDs that should be deleted.
                TR = triangulation(t,p) ;
                elock = edgeAttachments(TR,egfix) ;
                tq = gettrimeshquan(p,t);
                elock = unique(cell2mat(elock'));
                dmy = elock(tq.qm(elock) < 0.25);
                badtria = t(dmy,:);
                del     = badtria(badtria > nfix) ;
            end

        end % end mesh generator


    end % end methods


end % end class