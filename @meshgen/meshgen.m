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
    %
    %   Available options:
    %         ef            % edgefx class
    %         bou           % geodata class
    %         h0            % minimum edge length (optional if bou exists)
    %         bbox          % bounding box [xmin,ymin; xmax,ymax] (manual specification, no bou)
    %         proj          % structure containing the m_map projection info
    %         plot_on       % flag to plot (def: 1) or not (0)
    %         nscreen       % how many it to plot and write temp files (def: 5)
    %         itmax         % maximum number of iterations.
    %         pfix          % fixed node positions (nfix x 2 )
    %         egfix         % edge constraints
    %         outer         % meshing boundary (manual specification, no bou)
    %         inner         % island boundaries (manual specification, no bou)
    %         mainland      % the shoreline boundary (manual specification, no bou)
    %         fixboxes      % a flag that indicates which boxes will use fixed constraints
    %         memory_gb     % memory in GB allowed to use for initial rejector
    %         cleanup       % logical flag or string to trigger cleaning of topology (default is on).
    %         direc_smooth  % logical flag to trigger direct smoothing of mesh in the cleanup
    %         dj_cutoff     % the cutoff area fraction for disjoint portions to delete
    %         qual_tol      % tolerance for the accepted negligible change in quality
    %         enforceWeirs  % whether or not to enforce weirs in meshgen
    %         enforceMin    % whether or not to enfore minimum edgelength for all edgefxs
    % delaunay_elim_on_exit % whether or not to run delaunay_elim on exit of meshgen
    % improve_with_reduced_quality % whether or not to allow mesh improvements with decreases in mesh quality
    %      improve_boundary % run a gradient desc. to improve the boundary conformity
    %       high_fidelity   % flag to form pfix and egfix for this domain
    %
    properties
        fd            % handle to distance function
        fh            % handle to edge function
        h0            % minimum edge length
        edgefx        % edgefx class
        bbox          % bounding box [xmin,ymin; xmax,ymax]
        pfix          % fixed node positions (nfix x 2 )
        egfix         % edge constraints
        fixboxes      % a flag that indicates which boxes will use fixed constraints
        plot_on       % flag to plot (def: 1) or not (0)
        nscreen       % how many it to plot and write temp files (def: 5)
        bou           % geodata class
        ef            % edgefx class
        itmax         % maximum number of iterations.
        outer         % meshing boundary (manual specification, no bou)
        inner         % island boundaries (manual specification, no bou)
        mainland      % the shoreline boundary (manual specification, no bou)
        boubox        % the bbox as a polygon 2-tuple
        inpoly_flip   % used to flip the inpoly test to determine the signed distance.
        memory_gb     % memory in GB allowed to use for initial rejector
        cleanup       % logical flag or string to trigger cleaning of topology (default is on).
        direc_smooth  % logical flag to trigger direct smoothing of mesh in the cleanup
        dj_cutoff     % the cutoff area fraction for disjoint portions to delete
        grd = msh();  % create empty mesh class to return p and t in.
        qual          % mean, lower 3rd sigma, and the minimum element quality.
        qual_tol      % tolerance for the accepted negligible change in quality
        proj          % structure containing the m_map projection info
        anno          % Approx. Nearest Neighbor search object.
        annData       % datat contained with KD-tree in anno
        Fb            % bathymetry data interpolant
        enforceWeirs  % whether or not to enforce weirs in meshgen
        enforceMin    % whether or not to enfore minimum edgelength for all edgefxs
        improve_boundary  % improve the boundary representation
        high_fidelity     % flag to form pfix and egfix for this domain
        delaunay_elim_on_exit % whether or not to run delaunay_elim on exit of meshgen
        improve_with_reduced_quality % whether or not to allow mesh improvements with decreases in mesh quality
    end
    
    
    
    
    methods
        
        
        function obj = plot(obj)
            if ~isempty(obj.pfix) && ~isempty(obj.egfix)
                figure; hold on; % Initialize figure and hold for multiple plots
                
                % Plot the bounding boxes for visual aid
                for box_number = 1:length(obj.boubox)
                    iboubox = obj.boubox{box_number};
                    plot(iboubox(:,1),iboubox(:,2), 'g-', 'LineWidth', 2, 'DisplayName', 'Bounding Box');
                    
                    % Plot outer boundaries
                    touter = obj.outer(box_number);
                    if ~isempty(touter)
                        tedges = Get_poly_edges(touter{1}); % Assuming Get_poly_edges is defined elsewhere
                        [touter, ~] = filter_polygon_constraints(touter{1}, tedges, obj.boubox, box_number);
                        plot(touter(:,1), touter(:,2), 'bx', 'DisplayName', 'Outer Boundary');
                    end
                    
                    % Plot inner boundaries
                    tinner = obj.inner(box_number);
                    if ~isempty(tinner) && ~isempty(tinner{1})
                        tedges = Get_poly_edges(tinner{1}); % Assuming Get_poly_edges is defined elsewhere
                        [tinner, ~] = filter_polygon_constraints(tinner{1}, tedges, obj.boubox, box_number);
                        plot(tinner(:,1), tinner(:,2), 'rx', 'DisplayName', 'Inner Boundary');
                    end
                end
                
                % Plot the fixed constraints as squares and edges in black
                if exist('drawedge2', 'file') == 2 % Check if drawedge2 exists
                    drawedge2(obj.pfix, obj.egfix, [0,0,0]); % Assuming drawedge2 is a custom function
                else
                    % Alternative plotting if drawedge2 is not available
                    plot(obj.pfix(:,1), obj.pfix(:,2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Fixed Points');
                    for i = 1:size(obj.egfix, 1)
                        plot(obj.pfix(obj.egfix(i,:), 1), obj.pfix(obj.egfix(i,:), 2), 'k-', 'LineWidth', 2, 'DisplayName', 'Fixed Edges');
                    end
                end
                
                axis equal;
                title('Mesh Generation Constraints are black lines');
                xlabel('Longitude');
                ylabel('Latitude');
                legend('show', 'Location', 'bestoutside');
            else
                disp('No constraints to plot!');
            end
        end
        
        
        
        
        % class constructor/default grd generation options
        function obj = meshgen(varargin)
            % Check for m_map dir
            M_MAP_EXISTS=0;
            if exist('m_proj','file')==2
                M_MAP_EXISTS=1 ;
            end
            if M_MAP_EXISTS~=1
                error('Where''s m_map? Chief, you need to read the user guide')
            end
            
            
            % Check for utilties dir
            UTIL_DIR_EXISTS=0 ;
            if exist('inpoly.m','file')
                UTIL_DIR_EXISTS=1;
            end
            if UTIL_DIR_EXISTS~=1
                error('Where''s the utilities directory? Chief, you need to read the user guide')
            end
            
            
            p = inputParser;
            % unpack options and set default ones, catch errors.
            
            
            defval = 0; % placeholder value if arg is not passed.
            % add name/value pairs
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
            addOptional(p,'enforceWeirs',1);
            addOptional(p,'enforceMin',1);
            addOptional(p,'delaunay_elim_on_exit',1);
            addOptional(p,'improve_with_reduced_quality',0);
            
            
            % parse the inputs
            parse(p,varargin{:});
            
            
            %if isempty(varargin); return; end
            % store the inputs as a struct
            inp = p.Results;
     
            
            % kjr...order these argument so they are processed in a predictable
            % manner. Process the general opts first, then the OceanMesh
            % classes...then basic non-critical options.
            inp = orderfields(inp,{'h0','bbox','enforceWeirs','enforceMin',...
                'delaunay_elim_on_exit','improve_with_reduced_quality',...
                'fh',...
                'inner','outer','mainland',...
                'bou','ef',... %<--OceanMesh classes come after
                'egfix','pfix','fixboxes',...
                'plot_on','nscreen','itmax',...
                'memory_gb','qual_tol','cleanup',...
                'direc_smooth','dj_cutoff',...
                'big_mesh','proj'});
            % get the fieldnames of the edge functions
            fields = fieldnames(inp);
            % loop through and determine which args were passed.
            % also, assign reasonable default values if some options were
            % not assigned.
            for i = 1 : numel(fields)
                type = fields{i};
                switch type
                    % parse aux options first
                    case('h0')
                        % Provide in meters
                        obj.h0 = inp.(fields{i});
                    case('fh')
                        if isa(inp.(fields{i}),'function_handle')
                            obj.fh = inp.(fields{i});
                        end
                        % can't check for errors here yet.
                    case('bbox')
                        obj.bbox= inp.(fields{i});
                        if iscell(obj.bbox)
                            % checking bbox extents
                            ob_min = obj.bbox{1}(:,1);
                            ob_max = obj.bbox{1}(:,2);
                            for ii = 2:length(obj.bbox)
                                if any(obj.bbox{ii}(:,1) < ob_min) || ...
                                        any(obj.bbox{ii}(:,2) > ob_max)
                                    error(['Outer bbox must contain all ' ...
                                        'inner bboxes: inner box #' ...
                                        num2str(ii) ' violates this'])
                                end
                            end
                        end
                        
                        
                        % if user didn't pass anything explicitly for
                        % bounding box make it empty so it can be populated
                        % from ef as a cell-array
                        if obj.bbox(1)==0
                            obj.bbox = [];
                        end
                    case('pfix')
                        obj.pfix= inp.(fields{i});
                        if obj.pfix(1)~=0
                            obj.pfix(:,:) = inp.(fields{i});
                        else
                            obj.pfix = [];
                        end
                        if obj.enforceWeirs
                            for j = 1 : length(obj.bou)
                                if  ~isempty(obj.bou{j}.weirPfix)
                                    obj.pfix = [obj.pfix ; obj.bou{j}.weirPfix];
                                end
                            end
                        end
                    case('egfix')
                        obj.egfix= inp.(fields{i});
                        if ~isempty(obj.egfix) && obj.egfix(1)~=0
                            obj.egfix = inp.(fields{i});
                        else
                            obj.egfix = [];
                        end
                        if obj.enforceWeirs
                            for j = 1 : length(obj.bou)
                                if ~isempty(obj.bou{j}.weirEgfix) && ~isempty(obj.egfix)
                                    obj.egfix = [obj.egfix ; obj.bou{j}.weirEgfix+max(obj.egfix(:))];
                                elseif isempty(obj.egfix)
                                    obj.egfix =  obj.bou{j}.weirEgfix;
                                end
                            end
                        end
                        obj.egfix = renumberEdges(obj.egfix);
                    case('fixboxes')
                        obj.fixboxes= inp.(fields{i});
                    case('bou')
                        % got it from user arg
                        if obj.outer~=0, continue; end
                        
                        obj.outer = {} ;
                        obj.inner = {} ;
                        obj.mainland = {} ;
                        
                        obj.bou = inp.(fields{i});
                        
                        % handle when not a cell
                        if ~iscell(obj.bou)
                            boutemp = obj.bou;
                            obj.bou = cell(1);
                            obj.bou{1} = boutemp;
                        end
                        
                        % then the geodata class was provide, unpack
                        for ee = 1:length(obj.bou)
                            try
                                arg = obj.bou{ee} ;
                            catch
                                arg = obj.bou;
                            end
                            if isa(arg,'geodata')
                                
                                obj.high_fidelity{ee} = obj.bou{ee}.high_fidelity;
                                                                
                                obj.outer{ee} = obj.bou{ee}.outer;
                                obj.inner{ee} = obj.bou{ee}.inner;
                                
                                % save bathy interpolant to meshgen
                                if ~isempty(obj.bou{ee}.Fb)
                                    obj.Fb{ee} = obj.bou{ee}.Fb ;
                                end
                                
                                if ~isempty(obj.inner{ee}) && ...
                                        obj.inner{ee}(1)~= 0
                                    obj.outer{ee} = [obj.outer{ee};
                                        obj.inner{ee}];
                                end
                                obj.mainland{ee} = obj.bou{ee}.mainland;
                                obj.boubox{ee} = obj.bou{ee}.boubox;
                                obj.inpoly_flip{ee} = obj.bou{ee}.inpoly_flip;

                            end
                        end
                        
                        
                    case('ef')
                        tmp = inp.(fields{i});
                        if isa(tmp, 'function_handle')
                            error('Please specify your edge function handle through the name/value pair fh');
                        end
                        obj.ef = tmp;
                        
                        
                        % handle when not a cell
                        if ~iscell(obj.ef)
                            eftemp = obj.ef;
                            obj.ef = cell(1);
                            obj.ef{1} = eftemp;
                        end
                        
                        
                        % Gather boxes from ef class.
                        for ee = 1 : length(obj.ef)
                            if isa(obj.ef{ee},'edgefx')
                                obj.bbox{ee} = obj.ef{ee}.bbox;
                            end
                        end
                        
                        
                        % checking bbox extents
                        if iscell(obj.bbox)
                            ob_min = obj.bbox{1}(:,1);
                            ob_max = obj.bbox{1}(:,2);
                            for ii = 2:length(obj.bbox)
                                if any(obj.bbox{ii}(:,1) < ob_min) || ...
                                        any(obj.bbox{ii}(:,2) > ob_max)
                                    error(['Outer bbox must contain all ' ...
                                        'inner bboxes: inner box #' ...
                                        num2str(ii) ' violates this'])
                                end
                            end
                        end
                        
                        
                        % kjr 2018 June: get h0 from edge functions
                        for ee = 1:length(obj.ef)
                            if isa(obj.ef{ee},'edgefx')
                                obj.h0(ee) = obj.ef{ee}.h0;
                            end
                        end
                        
                        
                        % kjr 2018 smooth the outer automatically
                        if length(obj.ef) > 1
                            % kjr 2020, ensure the min. sizing func is
                            % used
                            if obj.enforceMin
                                obj.ef = enforce_min_ef(obj.ef);
                            end
                            obj.ef = smooth_outer(obj.ef,obj.Fb);
                        end
                        
                        
                        % Save the ef interpolants into the edgefx
                        for ee = 1:length(obj.ef)
                            if isa(obj.ef{ee},'edgefx')
                                obj.fh{ee} = @(p)obj.ef{ee}.F(p);
                            end
                        end
                        
                        
                    case('plot_on')
                        obj.plot_on= inp.(fields{i});
                    case('nscreen')
                        obj.nscreen= inp.(fields{i});
                        if obj.nscreen ~=0
                            obj.nscreen = inp.(fields{i});
                            obj.plot_on = 1;
                        else
                            obj.nscreen = 5; % default
                        end
                    case('itmax')
                        obj.itmax= inp.(fields{i});
                        if obj.itmax ~=0
                            obj.itmax = inp.(fields{i});
                        else
                            obj.itmax = 100;
                            warning('No itmax specified, itmax set to 100');
                        end
                    case('qual_tol')
                        obj.qual_tol = inp.(fields{i});
                        if obj.qual_tol ~=0
                            obj.qual_tol = inp.(fields{i});
                        else
                            obj.qual_tol = 0.01;
                        end
                    case('inner')
                        if ~isa(obj.bou,'geodata')
                            obj.inner = inp.(fields{i});
                        end
                    case('outer')
                        if ~isa(obj.bou,'geodata')
                            obj.outer = inp.(fields{i});
                            if obj.inner(1)~=0
                                obj.outer = [obj.outer; obj.inner];
                            end
                        end
                    case('mainland')
                        if ~isa(obj.bou,'geodata')
                            obj.mainland = inp.(fields{i});
                        end
                    case('memory_gb')
                        if ~isa(obj.bou,'memory_gb')
                            obj.memory_gb = inp.(fields{i});
                        end
                    case('cleanup')
                        obj.cleanup = inp.(fields{i});
                        if isempty(obj.cleanup) || obj.cleanup == 0
                            obj.cleanup = 'none';
                        elseif obj.cleanup == 1
                            obj.cleanup = 'default';
                        end
                    case('dj_cutoff')
                        obj.dj_cutoff = inp.(fields{i});
                    case('direc_smooth')
                        obj.direc_smooth = inp.(fields{i});
                    case('proj')
                        obj.proj = inp.(fields{i});
                        % default CPP
                        if obj.proj == 0; obj.proj = 'equi'; end
                        if ~isempty(obj.bbox)
                            % kjr Oct 2018 use outer coarsest box for
                            % multiscale meshing
                            lon_mi = obj.bbox{1}(1,1)-obj.h0(1)/1110;
                            lon_ma = obj.bbox{1}(1,2)+obj.h0(1)/1110;
                            lat_mi = obj.bbox{1}(2,1)-obj.h0(1)/1110;
                            lat_ma = obj.bbox{1}(2,2)+obj.h0(1)/1110;
                        else
                            lon_mi = -180; lon_ma = 180;
                            lat_mi = -90; lat_ma = 90;
                        end
                        % Set up projected space
                        dmy = msh() ;
                        dmy.p(:,1) = [lon_mi; lon_ma];
                        dmy.p(:,2) = [lat_mi; lat_ma];
                        del = setProj(dmy,1,obj.proj) ;
                    case('enforceWeirs')
                        obj.enforceWeirs = inp.(fields{i});
                    case('enforceMin')
                        obj.enforceMin = inp.(fields{i});
                    case('delaunay_elim_on_exit')
                        obj.delaunay_elim_on_exit = inp.(fields{i});
                    case('improve_with_reduced_quality')
                        obj.improve_with_reduced_quality = inp.(fields{i});
                end
            end
            
            
            if isempty(varargin); return; end
            
            
            % error checking
            if isempty(obj.boubox) && ~isempty(obj.bbox)
                % Make the bounding box 5 x 2 matrix in clockwise order if
                % it isn't present. This case must be when the user is
                % manually specifying the PSLG.
                obj.boubox{1} = [obj.bbox(1,1) obj.bbox(2,1);
                    obj.bbox(1,1) obj.bbox(2,2); ...
                    obj.bbox(1,2) obj.bbox(2,2);
                    obj.bbox(1,2) obj.bbox(2,1); ...
                    obj.bbox(1,1) obj.bbox(2,1); NaN NaN];
            end
            if any(obj.h0==0), error('h0 was not correctly specified!'), end
            if isempty(obj.outer), error('no outer boundary specified!'), end
            if isempty(obj.bbox), error('no bounding box specified!'), end
            obj.fd = @dpoly;  % <-default distance fx accepts p and pv (outer polygon).
            % kjr build ANN object into meshgen
            obj = createANN(obj) ;
            
            
            global MAP_PROJECTION MAP_COORDS MAP_VAR_LIST
            obj.grd.proj    = MAP_PROJECTION ;
            obj.grd.coord   = MAP_COORDS ;
            obj.grd.mapvar  = MAP_VAR_LIST ;
            
            % Check and update geometry for high fidelity option
            tpfix = [];
            tegfix = [];
            for box_num = 1:length(obj.h0)

                
                % High fidelity - formation of point & edge constraints
                if obj.high_fidelity{box_num}

                    % Define polygons based on available geometry data
                    ml = obj.mainland{box_num};
                    il = obj.inner{box_num};
                    polys = {};
                    if ~isempty(il), polys{end+1} = il; end
                    if ~isempty(ml), polys{end+1} = ml; end
                    
                    if obj.cleanup == 1
                        warning('Setting cleanup to 0 since high_fidelity mode is on');
                        obj.cleanup = 0;
                    end
                    disp(['Redistributing vertices for box #' num2str(box_num)]);
                    
                    % Convert polygon data, removing NaNs and duplicates
                    poly = cell2mat(polys');
                    D = nandelim_to_cell(poly); % Custom function to handle NaN-delimited cells
                    D = D(~cellfun(@(p) all(isnan(p(:))), D)); % Remove cells with all NaNs
                    areas = cellfun(@(d) polyarea(d(:, 1), d(:, 2)), D);
                    [~, uniqueIdx] = uniquetol(areas, 1e-12); % Remove duplicate polygons
                    D = D(uniqueIdx);
                    
                    % Process each polygon
                    for i = 1:length(D)
                        my_poly = D{i};
                        if size(my_poly, 1) > 2
                            tedges = Get_line_edges(my_poly);
                            [pts, bnde] = filter_polygon_constraints(my_poly, tedges, obj.boubox, box_num);
                            if isempty(bnde), continue; end
                            poly_split = extdom_polygon(bnde, pts, 0, 1);
                            
                            for ii = 1:length(poly_split)
                                points = unique(poly_split{ii}, 'rows', 'stable');
                                % Option 2 simply constrains the vector as
                                % it is in the file.
                                if obj.high_fidelity{box_num} == 2
                                    tmp_pfix = points;
                                    tmp_egfix = Get_poly_edges([points; NaN,NaN]);
                                % Option 1: resamples based on the local
                                % mesh resolution
                                elseif obj.high_fidelity{box_num} == 1
                                    [tmp_pfix, tmp_egfix] = mesh1d(points, obj.fh, obj.h0./111e3,...
                                        [], obj.boubox, box_num, []);
                                end
                                
                                if size(tmp_pfix, 1) > 2
                                    [tmp_pfix, tmp_egfix] = fixgeo2(tmp_pfix, tmp_egfix);
                                    if max(tmp_egfix(:)) ~= size(tmp_pfix, 1), continue; end
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
                        
            % Final geometry check and update
            [tpfix, tegfix] = fixgeo2(tpfix, tegfix);
            checkfixp = setdiff(tpfix, fixmesh(tpfix), 'rows');
            if ~isempty(checkfixp), error('Duplicate fixed points detected, cannot proceed'); end
            
            % Update object with user-passed fixed points
            obj.pfix = [obj.pfix; tpfix];
            
            if isempty(obj.egfix)
                obj.egfix = tegfix;
            else
                obj.egfix = [obj.egfix; tegfix + max(obj.egfix(:))];
            end
            
            if ~isempty(obj.pfix)
                disp(['Using ', num2str(length(obj.pfix)), ' fixed points.']);
            end
            
            if ~isempty(obj.egfix)
                if max(obj.egfix(:)) > length(obj.pfix)
                    error('FATAL: egfix does not index correctly into pfix.');
                end
                disp(['Using ', num2str(length(obj.egfix)), ' fixed edges.']);
                warning('Please check fixed constraints: plot(mshopts)');
            end
        end
        
        
        % Creates Approximate nearest neighbor objects on start-up
        function  obj = createANN(obj)
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
                    % This line removes the line that can appear in the
                    % center for a global mesh
                    dataset(abs(dataset(:,1)) > 180-1e-6,:) = [];
                    dataset(abs(dataset(:,1)) < 1e-6,:) = [];
                end
                [dataset(:,1),dataset(:,2)] = m_ll2xy(dataset(:,1),dataset(:,2));
                dataset(isnan(dataset(:,1)),:) = [];
                dmy = ann(dataset');
                obj.anno{box_num} = dmy;
                obj.annData{box_num}=dataset;
            end
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

