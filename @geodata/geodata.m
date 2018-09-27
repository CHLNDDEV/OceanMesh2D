classdef geodata
    %   GEODATA: Geographical data class
    %   Handles geographical data describing coastlines or other features in
    %   the form of a shapefile and topobathy in the form of a DEM
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
        bbox  % bounding coordinates (same as mesh and grdgen classes)
        boubox % bbox coordinates as a polygon
        mainland % mainland boundary.
        outer % outer boundary.
        inner % islands.
        inpoly_flip
        contourfile %  cell-array of shapefile
        demfile     % filename of dem
        h0 % min. edgelength.
        window=5 % smoothing window on boundary
        fp
        Fb  % linear gridded interpolant of DEM
        x0y0 % bottom left of structure grid or position (0,0)
    end
    
    methods
        
        % class constructor/parse shapefile
        function obj = geodata(varargin)
            
            p = inputParser;
            defval = 0; % placeholder value if arg is not passed.
            % add name/value pairs
            addOptional(p,'bbox',defval);
            addOptional(p,'shp',defval);
            addOptional(p,'h0',defval);
            addOptional(p,'dem',defval);
            addOptional(p,'fp',defval);
            addOptional(p,'outer',defval);
            addOptional(p,'mainland',defval);
            
            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp=p.Results;
            % get the fieldnames of the edge functions
            fields = fieldnames(inp);
            % loop through and determine which args were passed.
            % also, assign reasonable default values if some options were
            % not assigned.
            for i = 1 : numel(fields)
                type = fields{i};
                switch type
                    % parse aux options first
                    case('bbox')
                        obj.bbox = inp.(fields{i});
                        if any(obj.bbox ~=0)
                            obj.bbox = inp.(fields{i});
                        elseif ~ischar(inp.(fields{contains(fields,'dem')}))
                            error('No bbox specified!');
                        end
                    case('h0')
                        obj.h0= inp.(fields{i});
                        if obj.h0 ~=0
                            obj.h0 = inp.(fields{i});
                        else
                            error('No h0 specified!');
                        end
                    case('shp')
                        obj.contourfile = inp.(fields{i});
                        if ~iscell(obj.contourfile) && ~ischar(obj.contourfile)
                            obj.contourfile = [];
                        end
                    case('dem')
                        obj.demfile= inp.(fields{i});
                        if obj.demfile ~=0
                            obj.demfile = inp.(fields{i});
                        else
                            obj.demfile = [];
                        end
                    case('fp')
                        obj.fp= inp.(fields{i});
                        if obj.fp ~=0
                            obj.fp = inp.(fields{i});
                        end
                    case('outer')
                        obj.outer = inp.(fields{i});
                        if obj.outer(1) ~=0
                            obj.outer = inp.(fields{i});
                        end
                    case('mainland')
                        obj.mainland = inp.(fields{i});
                        if obj.mainland(1) ~=0
                            obj.mainland = inp.(fields{i});
                        end
                end
            end
            
            if size(obj.bbox,1) == 1
                % get bbox from demfile
                try
                   x = double(ncread(obj.demfile,'lon'));
                   y = double(ncread(obj.demfile,'lat'));
                catch
                   x = double(ncread(obj.demfile,'x'));
                   y = double(ncread(obj.demfile,'y'));    
                end
                obj.bbox = [min(x) max(x); min(y) max(y)];
            end
            
            if ~isempty(obj.contourfile)
                centroid     = mean(obj.bbox(2,:));
                gridspace =    abs(obj.h0)/(cosd(centroid)*111e3);
                % Read polygon from shape file and make sure spacing is h0.;
                if ~iscell(obj.contourfile)
                    obj.contourfile = {obj.contourfile};
                end
                polygon_struct = Read_shapefile( obj.contourfile, [], obj.bbox, ...
                    gridspace, 0 );
                obj.mainland = polygon_struct.mainland;
                obj.outer    = polygon_struct.outer;
                obj.inner    = polygon_struct.inner;
                
                % for mainland
                if ~isempty(obj.mainland)
                    [la,lo]=my_interpm(obj.mainland(:,2),obj.mainland(:,1),gridspace/2);
                    obj.mainland = [];
                    obj.mainland(:,1) = lo; obj.mainland(:,2) = la;
                end
                
                % make sure it has equal spacing of h0/2
                [la,lo]=my_interpm(obj.outer(:,2),obj.outer(:,1),gridspace/2);
                obj.outer = [];
                obj.outer(:,1) = lo; obj.outer(:,2) = la;
                
                % for islands
                if ~isempty(obj.inner)
                    [la,lo]=my_interpm(obj.inner(:,2),obj.inner(:,1),gridspace/2);
                    obj.inner = [];
                    obj.inner(:,1) = lo; obj.inner(:,2) = la;
                end
                clearvars lo la
                
                % smooth coastline
                if ~isempty(obj.outer)
                    obj.outer = smooth_coastline(obj.outer,obj.window,0);
                end
                if ~isempty(obj.mainland)
                    obj.mainland = smooth_coastline(obj.mainland,obj.window,0);
                end
                if ~isempty(obj.inner)
                    obj.inner = smooth_coastline(obj.inner,obj.window,0);
                end
                
                % WJP: Jan 25, 2018 check the polygon
                obj = check_connectedness_inpoly(obj);
                
                % Make the bounding box 5 x 2 matrix in clockwise order
                obj.boubox = [obj.bbox(1,1) obj.bbox(2,1);
                    obj.bbox(1,1) obj.bbox(2,2); ...
                    obj.bbox(1,2) obj.bbox(2,2);
                    obj.bbox(1,2) obj.bbox(2,1); ...
                    obj.bbox(1,1) obj.bbox(2,1); NaN NaN];
                iboubox = obj.boubox;
                iboubox(:,1) = 1.10*iboubox(:,1)+(1-1.10)*mean(iboubox(1:end-1,1));
                iboubox(:,2) = 1.10*iboubox(:,2)+(1-1.10)*mean(iboubox(1:end-1,2));
                
                % KJR: May 13, 2018 coarsen portions of outer, mainland 
                % and inner outside bbox (made into function WP)
                obj.outer = coarsen_polygon(obj.outer,iboubox);
              
                if ~isempty(obj.inner)
                    obj.inner = coarsen_polygon(obj.inner,iboubox);
                end
                
                if ~isempty(obj.mainland)
                    obj.mainland = coarsen_polygon(obj.mainland,iboubox);
                end
                
                disp(['Read in contourfile ',obj.contourfile]);
                
            else
                % then user defined file passed
                centroid     = mean(obj.bbox(2,:));
                gridspace =    abs(obj.h0)/(cosd(centroid)*111e3);
                % make sure it has equal spacing of h0/2
                [la,lo]=my_interpm(obj.outer(:,2),obj.outer(:,1),gridspace/2);
                obj.outer = [];
                obj.outer(:,1) = lo; obj.outer(:,2) = la;
                if isempty(obj.outer)
                    error('Outer segment is required to mesh!');
                end
                % for mainland
                if ~isempty(obj.mainland)
                    [la,lo]=my_interpm(obj.mainland(:,2),obj.mainland(:,1),gridspace/2);
                    obj.mainland = [];
                    obj.mainland(:,1) = lo; obj.mainland(:,2) = la;
                else
                    error('Mainland segment is required to mesh!');
                end
                % for islands
                if ~isempty(obj.inner)
                    [la,lo]=my_interpm(obj.inner(:,2),obj.inner(:,1),gridspace/2);
                    obj.inner = [];
                    obj.inner(:,1) = lo; obj.inner(:,2) = la;
                end
                clearvars lo la
                
                % Make the bounding box 5 x 2 matrix in clockwise order
                obj.boubox = [obj.bbox(1,1) obj.bbox(2,1);
                    obj.bbox(1,1) obj.bbox(2,2); ...
                    obj.bbox(1,2) obj.bbox(2,2);
                    obj.bbox(1,2) obj.bbox(2,1); ...
                    obj.bbox(1,1) obj.bbox(2,1); NaN NaN];
                
            end
            
            if ~isempty(obj.demfile)
                try
                   x = double(ncread(obj.demfile,'lon'));
                   y = double(ncread(obj.demfile,'lat'));
                catch
                   x = double(ncread(obj.demfile,'x'));
                   y = double(ncread(obj.demfile,'y'));    
                end
                % Read with Buffer
                centroid  = mean(obj.bbox(2,:));
                gridspace = abs(obj.h0)/(cosd(centroid)*111e3);
                I = find(x >= obj.bbox(1,1) & ...
                         x <= obj.bbox(1,2));
                J = find(y >= obj.bbox(2,1) & ...
                         y <= obj.bbox(2,2));
                x = x(I); y = y(J);
                % Find name of z value (use first one that has 2 dimensions
                finfo = ncinfo(obj.demfile);
                for ii = 1:length(finfo.Variables)
                    if length(finfo.Variables(ii).Size) == 2
                        zvarname = finfo.Variables(ii).Name;
                        break
                    end
                end
                demz = single(ncread(obj.demfile,zvarname,...
                              [I(1) J(1)],[length(I) length(J)]));
                obj.Fb   = griddedInterpolant({x,y},demz,...
                    'linear','nearest');
                obj.x0y0 = [x(1),y(1)];
                clear x y demz
                disp(['Read in demfile ',obj.demfile]);
            end
            
            % then no dem
            if isempty(obj.x0y0)
                obj.x0y0 = [obj.bbox(1,1), obj.bbox(2,1)];
            end
            
            function polyobj = coarsen_polygon(polyobj,iboubox)
                % coarsen a polygon object
                idx = find(isnan(polyobj(:,1)));
                C = mat2cell(polyobj,diff([0;idx]),2);
                new = cell(length(C),1);
                % for each segment
                for iii = 1 : length(C)
                    j = 0 ; k = 0;
                    [in] = inpoly(C{iii},iboubox(1:end-1,:));
                    % Initialise for speed
                    new{iii} = zeros(length(C{iii}),2);
                    while j < length(C{iii})-1
                        j = j + 1 ;
                        if in(j) % point is in domain keep it
                            k = k + 1;
                            new{iii}(k,:) = C{iii}(j,:);
                        else % pt is out of domain
                            bd = min([j+200,length(in)-1]) ;
                            exte = min(200,bd - j);
                            if sum(in(j:bd))==0 % if next hundred points are out, then we can decimate
                                k = k + 1 ;
                                new{iii}(k,:) = C{iii}(j,:);
                                k = k + 1 ;
                                new{iii}(k,:) = C{iii}(j+exte,:);
                                j = j + exte ;
                            else % otherwise keep
                                k = k + 1 ;
                                new{iii}(k,:) = C{iii}(j,:);
                            end

                        end
                    end
                    new{iii}(k+1,:) = [NaN NaN] ;
                    new{iii}(k+2:end,:) = [];
                end
                polyobj = cell2mat(new);
            end
            
        end
        
        % WJP: Check for connected polygons and if not connected,
        %      whether to flip the inpoly result to ensure meshing
        %      on the coastal side of the polygon
        function obj = check_connectedness_inpoly(obj)
            centroid     = mean(obj.bbox(2,:));
            gridspace =    abs(obj.h0)/(cosd(centroid)*111e3);
            
            obj.inpoly_flip = 0;
            % return if outer polygon is connected
            shpEnd = find(isnan(obj.outer(:,1)));
            [~,loc] = max(diff(shpEnd));
            shpEnd = vertcat(0,shpEnd); loc = loc+1;
            if abs(sum(obj.outer(shpEnd(loc)+1,:)) - ...
                    sum(obj.outer(shpEnd(loc+1)-1,:))) < eps
                return;
            end
            disp('Warning: Shapefile is unconnected... continuing anyway')
            % outer polygon is not connected, check for inpoly goodness
            
            % read the GSHHS checker
            ps = Read_shapefile( {'GSHHS_l_L1'}, [], ...
                obj.bbox, gridspace, 0 );
            
            % make a fake tester grid
            x = linspace(obj.bbox(1,1),obj.bbox(1,2),100);
            y = linspace(obj.bbox(2,1),obj.bbox(2,2),100);
            edges = Get_poly_edges( [ps.outer; ps.inner] );
            in_Test = inpoly([x',y'],[ps.outer; ps.inner],edges);
            edges = Get_poly_edges( [obj.outer; obj.inner] );
            in_Shpf = inpoly([x',y'],[obj.outer; obj.inner],edges);
            % if more than half of thepoints disagree between Test and Shpf
            % lets flip the inpoly
            if length(find(xor(in_Shpf,in_Test))) > 50
                obj.inpoly_flip = 1;
                disp(['Unconnected shapefile is inconsistent ' ...
                    'with GHSSS test file, flipping the inpoly test'])
            end
            
            % if flooplaind meshing, flip the inpoly test
            if obj.fp
                obj.inpoly_flip = mod(1,obj.inpoly_flip);
            end
        end

        % close geometric countour vectors by clipping with a boubox 
        % updates gdat.outer so that the meshing domain is correctly defined 
        function obj = close(obj,seed)
            % Clips the mainland segment with the boubox. 
            % Performs a breadth-first search given a seed position
            % of the piecewise-straight line graph (PSLG) that is used to define the meshing boundary.
            % This returns back an updated geodata class instance with the outer boundary clipped with the boubox.
            % kjr,und,chl 2018
            
            if(nargin < 2),error('Must supply coordinate seed to perform flood-fill'); end
            
            geom = [obj.mainland; obj.inner; obj.boubox] ;
            
            [NODE,PSLG]=getnan2(geom) ;
            
            [NODE2,PSLG2,PART2] = bfsgeo2(NODE,PSLG,seed) ;
            
            POLY = extdom_polygon(PSLG2(PART2{1},:),NODE2,-1) ;
            
            new_outer = cell2mat(POLY') ;
            
            obj.outer = new_outer ;
            
            % reset this to default
            obj.inpoly_flip = 0 ;
        end
        
        
        % plot shp object on projected map
        function plot(obj,type,projection)
            if nargin == 1
                type = 'shp';
            end
            if nargin < 3
                projection = 'Transverse Mercator';
            end
            bufx = 0.2*(obj.bbox(1,2) - obj.bbox(1,1));
            bufy = 0.2*(obj.bbox(2,2) - obj.bbox(2,1));
            m_proj(projection,...
                   'long',[obj.bbox(1,1) - bufx, obj.bbox(1,2) + bufx],...
                   'lat',[obj.bbox(2,1) - bufy, obj.bbox(2,2) + bufy]);
            switch type
                case('dem')
                    % interpolate DEM's bathy linearly onto our edgefunction grid.
                    [demx,demy] = ndgrid(obj.x0y0(1):obj.h0/111e3:obj.bbox(1,2), ...
                                     obj.x0y0(2):obj.h0/111e3:obj.bbox(2,2));
                    demz = obj.Fb(demx,demy);                     
                    edges = Get_poly_edges( [obj.outer; obj.inner] );
                    in = inpoly([demx(:),demy(:)],[obj.outer; obj.inner], edges);
                    cmap = cmocean('deep');
                    hold on; m_fastscatter(demx(in),demy(in),demz(in),cmap); 
                    cb = colorbar; ylabel(cb,'topo-bathy depth [m]')
                case('omega')
                    % hatch the meshing domain, Omega
                    [demx,demy] = ndgrid(obj.x0y0(1):obj.h0/111e3:obj.bbox(1,2), ...
                                     obj.x0y0(2):obj.h0/111e3:obj.bbox(2,2));
                    edges = Get_poly_edges( [obj.outer; obj.inner] ); 
                    in = inpoly([demx(:),demy(:)],[obj.outer; obj.inner], edges);   
                    long = demx(~in); lati = demy(~in);
                    hold on; m_hatch(obj.boubox(1:end-1,1),...
                    obj.boubox(1:end-1,2),'cross',45,0.05);
                    m_plot(long,lati,'.','Color','white')
            end
            if ~isempty(obj.mainland)
                h1 = m_plot(obj.mainland(:,1),obj.mainland(:,2),...
                    'r-','linewi',1); hold on;
            end
            if ~isempty(obj.inner)
                h2 = m_plot(obj.inner(:,1),obj.inner(:,2),...
                    'g-','linewi',1); hold on;
            end
            m_plot(obj.boubox(:,1),obj.boubox(:,2),'k--','linewi',2,'MarkerSize',15);
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',10);
            legend([h1 h2],{'mainland' 'inner'},'Location','NorthWest')
            %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        end
        
    end
end

