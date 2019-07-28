classdef msh
    %   MSH: Mesh class
    %   Contains, handles, and builds properties of a mesh such as vertices,
    %   an element table, bathymetry and ADCIRC style boundary types
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
        title % mesh title
        p  % vertices
        t  % triangles
        b  % bathymetry
        bd % land boundaries
        op % open boundaries
        bx % slope of bathy in x direction
        by % slope of bathy in y direction
        f11 % A structure for the fort11 (initial density) values
        f13 % A struct of the fort13 attributes
        f15 % A struct of the fort15 inputs
        f19 % A struct for the fort19 non-periodic elevation bc
        f20 % A struct for the fort20 non-periodic flux/ele radiation bc
        f24 % A struct of the fort24 SAL values
        f2001 % A struct for the fort2001 non-periodic flux/ele sponge bc
        f5354 % A struct for the fort53001/54001 tidal ele/flux sponge bc
        proj   % Description of projected space (m_mapv1.4)
        coord  % Description of projected space (m_mapv1.4)
        mapvar % Description of projected space (m_mapv1.4)
        pfix   % fixed points that were constrained in msh
        egfix  % fixed eges that were constraind in msh
    end
    
    methods
        % constructor/read mesh into class.
        function obj = msh(fname,type)
            % Check for m_map dir
            M_MAP_EXISTS=0 ;
            if exist('m_proj','file')==2
              M_MAP_EXISTS=1 ;
            end
            if M_MAP_EXISTS~=1 
              error('Where''s m_map? Chief, you need to read the user guide')
            end

            % Check for utilties dir
            UTIL_DIR_EXISTS=0 ;
            if exist('inpoly.m','file')
              UTIL_DIR_EXISTS=1 ;
            end
            if UTIL_DIR_EXISTS~=1 
              error('Where''s the utilities directory? Chief, you need to read the user guide')
            end
            
            % Check for dataset dir
            DATASET_DIR_EXISTS=0 ;
            if exist('datasets','dir')==7
                DATASET_DIR_EXISTS=1 ;
            end
            if DATASET_DIR_EXISTS~=1
                warning('We suggest you to place your files in a datasets directory. Chief, you need to read the user guide')
            end

            if nargin == 0
                obj.title = 'OceanMesh2D';
                return
            end
            if nargin == 1
                type = '14';
            end
            if any(contains(type,'14'))
                bdflag = 1;
                if any(contains(type,'14nob'))
                    bdflag = 0;
                end
                [t,p,b,op,bd,title] = readfort14([fname '.14'],bdflag);
                obj.p  = p; obj.t  = t; obj.b  = b;
                obj.bd = bd; obj.op = op;
                obj.title = title;
            end
            if any(contains(type,'13'))
                obj.f13 = readfort13([fname '.13']);
            end
            if any(contains(type,'15'))
                if isempty(obj.op) ||  isempty(obj.bd)
                    error('Boundary data required to read f15...also read in f14.')
                end
                obj.f15 = readfort15([fname '.15'],obj.op,obj.bd);
            end
            if any(contains(type,'24'))
                if isempty(obj.p)
                    error('No vertices present to readfort24')
                end
                if isempty(obj.f15)
                    error('No f15 present to readfort24')
                end
                if obj.f15.ntif == 0
                    error('No constituents in f15 to readfort24')
                end
                obj.f24 = readfort24( [fname '.24'], obj.f15.ntif, ...
                    length(obj.p), {obj.f15.tipotag.name} );
            end
        end
        
        % write mesh to disk
        function write(obj,fname,type)
            if nargin == 1
                fname = 'fort_1';
            end
            if nargin < 3
                if isempty(obj.p)
                    error('No mesh, cannot write.')
                end
                
                % renumber it use RCM by default
                obj = renum(obj) ; 
                
                if isempty(obj.b)
                    b_t = 0*obj.p(:,1);
                else
                    b_t = obj.b;
                end
                writefort14( [fname '.14'], obj.t, obj.p, b_t, ...
                    obj.op , obj.bd ,obj.title ) ;
                if ~isempty(obj.f13)
                    writefort13( obj.f13, [fname '.13'] );
                end
                if ~isempty(obj.f15)
                    writefort15( obj.f15, [fname '.15'], obj.bd );
                end
                if ~isempty(obj.f24)
                    writefort24( obj.f24, [fname '.24'] );
                end
                if ~isempty(obj.f5354)
                    writefort5354( obj.f5354, fname );
                end
            else
                if any(contains(type,'14')) || any(contains(type,'ww3'))
                    if isempty(obj.p)
                        error('No mesh, cannot write.')
                    end
                    if isempty(obj.b)
                        b_t = 0*obj.p(:,1);
                    else
                        b_t = obj.b;
                    end
                    if any(contains(type,'14'))
                        writefort14( [fname '.14'] , obj.t, obj.p, b_t, ...
                            obj.op , obj.bd ,obj.title ) ;
                    end
                    if any(contains(type,'ww3'))
                        writeww3( [fname '.ww3'] , obj.t, obj.p, b_t, ...
                            obj.op , obj.title ) ;
                    end
                end
                if any(contains(type,'11')) && ~isempty(obj.f11)
                    writefort11( obj.f11, [fname '.11'] );
                end
                if any(contains(type,'13')) && ~isempty(obj.f13)
                    writefort13( obj.f13, [fname '.13'] );
                end
                if any(contains(type,'15')) && ~isempty(obj.f15)
                    writefort15( obj.f15, [fname '.15'], obj.bd );
                end
                if any(contains(type,'19')) && ~isempty(obj.f19)
                    writefort19( obj.f19, [fname '.19'] );
                end
                if any(contains(type,'2001')) && ~isempty(obj.f2001)
                    writefort19( obj.f2001, [fname '.2001'] );
                end
                if any(contains(type,'24')) && ~isempty(obj.f24)
                    writefort24( obj.f24, [fname '.24'] );
                end
                if any(contains(type,'5354')) && ~isempty(obj.f5354)
                    writefort5354( obj.f5354, fname );
                end
            end
        end
        
        % general plot function
        function h = plot(obj,type,proj,projtype,bou,numticks)
            % h = plot(obj,type,proj,projtype,bou,numticks)
            % 1) obj: msh object
            % 2) type: plot type, choose from:
            %    a) 'tri'  - (default) plots the triangulation
            %    b) 'bd'   - same as tri but with nodestrings plotted
            %    c) 'b'    - plots the bathymetry
            %    d) 'reso' - plots the element circumradius
            %    e) 'resodx' - plots the gradient in 'reso'
            %    f) 'slp'  - plots the bathymetric gradients
            %    g) 'itfric' - plots the internal_tide_friction values
            %    h) 'cfvals' - plots the quadratic bottom friction values
            %    additional --> 
            %    i)  add 'log' inside type to plot caxis in log space
            %    ii) add 'mesh' inside type to plot trimesh instead of trisurf
            % 3) proj: whether to plot in projected space or unprojected space
            %    a) 0       - plot in unprojected space
            %    b) 1       - plot in projected space (default)
            % 4) projtype: what projection to plot in if proj = 1. 
            %    default is the projection of the msh object
            % 5) bou: a local bounding box or polygon region to plot within
            % 6) numticks: number of tickmarks, if you plot in the 'log'
            %    space you can specify numticks = [numticks caxis_lower ...
            %                                      caxis_upper]
            if nargin < 2
                type = 'tri';
            end
            if nargin < 3
                proj = 1 ;
            end
            if nargin < 4
                projtype = [] ;
            end
            np_g = length(obj.p) ;
            if nargin < 6
                numticks = 10;
            end
            
            % kjr default behavior, just use what's in the .mat file
            if proj && isempty(projtype)
                global MAP_PROJECTION MAP_VAR_LIST MAP_COORDS
                if ~isempty(obj.coord)
                    % kjr 2018,10,17; Set up projected space imported from msh class
                    MAP_PROJECTION = obj.proj ;
                    MAP_VAR_LIST   = obj.mapvar ;
                    MAP_COORDS     = obj.coord ;
                    del = 0;
                    projtype = MAP_PROJECTION.name;
                end
            end
            
            % Handle user specified subdomain
            if nargin < 5 || isempty(bou)
                kept = (1:length(obj.p))';
            else
                if numel(bou) == 4
                    % i.e. is a bounding box
                    bou = [bou(1,1) bou(2,1);
                        bou(1,1) bou(2,2); ...
                        bou(1,2) bou(2,2);
                        bou(1,2) bou(2,1); ...
                        bou(1,1) bou(2,1)];
                end
                % Get a subset given by bou
                [obj,kept] = ExtractSubDomain(obj,bou);
            end
            
            % Set up projected space
            del = setProj(obj,proj,projtype) ;
            
            if del
                % This deletes any elements straddling the -180/180
                % boundary for plotting purposes
                xt = [obj.p(obj.t(:,1),1) obj.p(obj.t(:,2),1) ...
                     obj.p(obj.t(:,3),1) obj.p(obj.t(:,1),1)];
                dxt = diff(xt,[],2);
                obj.t(abs(dxt(:,1)) > 180 | abs(dxt(:,2)) > 180 | ...
                      abs(dxt(:,3)) > 180,:) = [];
            end
            
            logaxis = 0;
            idxl = strfind(type,'log');
            if ~isempty(idxl)
                logaxis = 1; type(idxl:idxl+2) = [];
            end
            mesh = 0;
            idxl = strfind(type,'mesh');
            if ~isempty(idxl)
                mesh = 1; type(idxl:idxl+3) = [];
            end    
            earthres = 0;
            if strcmp(type,'resoearth')
                type = 'reso'; earthres = 1;
            end
            
            switch type
                % parse aux options first
                case('tri')
                    figure;
                    if proj
                        m_triplot(obj.p(:,1),obj.p(:,2),obj.t);
                    else
                        simpplot(obj.p,obj.t);
                    end
                case('bd')
                    figure; hold on;
                    if proj
                        m_triplot(obj.p(:,1),obj.p(:,2),obj.t);
                    else
                        simpplot(obj.p,obj.t);
                    end
                    if ~isempty(obj.bd)
                        for nb = 1 : obj.bd.nbou
                            if obj.bd.ibtype(nb) == 94
                                if proj
                                    m_plot(obj.p(obj.bd.nbvv(:),1),...
                                        obj.p(obj.bd.nbvv(:),2),...
                                        'r.','linewi',1.2);
                                else
                                    plot(obj.p(obj.bd.nbvv(:),1),...
                                        obj.p(obj.bd.nbvv(:),2),...
                                        'r.','linewi',1.2);
                                end
                                % internal weirs
                            elseif obj.bd.ibtype(nb)  == 24
                                if proj
                                    % plot front facing
                                    m_plot(obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),2),'g-','linewi',1.2);
                                    % plot back facing
                                    m_plot(obj.p(obj.bd.ibconn(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.ibconn(1:obj.bd.nvell(nb),nb),2),'g-','linewi',1.2);
                                else
                                    % plot front facing
                                    plot(obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),2),'g-','linewi',1.2);
                                    % plot back facing
                                    plot(obj.p(obj.bd.ibconn(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.ibconn(1:obj.bd.nvell(nb),nb),2),'g','linewi',1.2);
                                end
                            elseif obj.bd.ibtype(nb)  == 20
                                if proj
                                    m_plot(obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),2),'r-','linewi',1.2);
                                else
                                    plot(obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),2),'r-','linewi',1.2);
                                end
                            else              
                                if proj
                                    m_plot(obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),2),'g-','linewi',1.2);
                                else
                                    plot(obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),1),...
                                        obj.p(obj.bd.nbvv(1:obj.bd.nvell(nb),nb),2),'g-','linewi',1.2);
                                end
                            end
                        end
                    end
                    if ~isempty(obj.op)
                        for nb = 1 : obj.op.nope
                            if proj
                                m_plot(obj.p(obj.op.nbdv(1:obj.op.nvdll(nb),nb),1),...
                                    obj.p(obj.op.nbdv(1:obj.op.nvdll(nb),nb),2),'b-','linewi',3.2);
                            else
                                plot(obj.p(obj.op.nbdv(1:obj.op.nvdll(nb),nb),1),...
                                    obj.p(obj.op.nbdv(1:obj.op.nvdll(nb),nb),2),'b-','linewi',1.2);
                            end
                        end
                    end
                    if isempty(obj.bd) && isempty(obj.op)
                        disp('bounbaries are empty!');
                    end
                case('b')
                    figure;
                    if logaxis
                        q = log10(max(1,abs(obj.b))); % plot on log scale
                    else
                        if exist('demcmap','file')
                            q = -obj.b;
                        else
                            q = obj.b;
                        end
                    end
                    if proj
                        if mesh
                            m_trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                        else
                            m_trisurf(obj.t,obj.p(:,1),obj.p(:,2),q);
                        end
                    else
                        if mesh
                            trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                        else
                            trisurf(obj.t,obj.p(:,1),obj.p(:,2),q)
                            shading interp;
                        end
                        view(2); 
                    end
                    if logaxis
                        cmocean('deep',numticks(1)-1); 
                    else
                        if exist('demcmap','file')
                            demcmap(q);
                        else   
                            cmocean('topo','pivot',min(max(q),0)); 
                        end
                    end
                    cb = colorbar;
                    if logaxis
                        if length(numticks) == 3
                            desiredTicks = round(10.^(linspace(...
                                           log10(numticks(2)),...
                                           log10(numticks(3)),...
                                           numticks(1))),1);
                        else
                            desiredTicks = round(10.^(linspace(min(q),...
                                                   max(q),numticks(1))),1);
                        end
                        caxis([log10(min(desiredTicks)) log10(max(desiredTicks))]);
                        cb.Ticks     = log10(desiredTicks);
                        for i = 1 : length(desiredTicks)
                            cb.TickLabels{i} = num2str(desiredTicks(i));
                        end
                    end
                    ylabel(cb,'m below geoid');
                    title('mesh topo-bathy');
                case('slp')
                    figure;
                    if proj
                        if mesh
                            m_trimesh(obj.t,obj.p(:,1),obj.p(:,2),...
                                      hypot(obj.bx,obj.by));  
                        else
                            m_trisurf(obj.t,obj.p(:,1),obj.p(:,2),...
                                      hypot(obj.bx,obj.by));
                        end
                    else
                        if mesh 
                            trimesh(obj.t,obj.p(:,1),obj.p(:,2),...
                                    hypot(obj.bx,obj.by));
                        else
                            trisurf(obj.t,obj.p(:,1),obj.p(:,2),...
                                    hypot(obj.bx,obj.by));
                            shading flat   
                        end
                        view(2);
                    end
                    colormap(cmocean('thermal'));
                    cb = colorbar; ylabel(cb,'slope');
                    caxis([0 0.25])
                case('ob') % outer boundary of mesh
                    [~,bpt] = extdom_edges2(obj.t,obj.p);
                    if proj
                        figure, m_plot(bpt(:,1),bpt(:,2),'r.');
                    else
                        figure, plot(bpt(:,1),bpt(:,2),'r.');
                    end
                case('reso')
                    % Get bar lengths
                    if earthres
                        [bars,barlen] = GetBarLengths(obj,0);
                        % sort bar lengths in ascending order
                        [barlen,IA] = sort(barlen,'descend');
                        bars = bars(IA,:);
                        % get the minimum bar length for each node
                        [B1,IB] = unique(bars(:,1),'last');
                        [B2,IC] = unique(bars(:,2),'last');
                        d1 = NaN*obj.p(:,1); d2 = NaN*obj.p(:,1);
                        d1(B1) = barlen(IB); d2(B2) = barlen(IC);  
                        z = min(d1,d2);
                    else
                        % get the points on the current projection
                        [X,Y]= m_ll2xy(obj.p(:,1),obj.p(:,2));
                        % get the circumcenter radius of each element
                        TR = triangulation(obj.t,X,Y);
                        [~,cr] = circumcenter(TR);
                        % Get the element connectivity
                        [vtoe,nne] = VertToEle(obj.t);
                        % Make sure null value adds zero contribution
                        cr(end+1) = 0;
                        vtoe(vtoe == 0) = length(obj.t) + 1;
                        % Sum up all element contributions to each node and
                        % divide by number of connected elements
                        z = sum(cr(vtoe))./nne;
                        % scale by earth radius
                        Re = 6378.137e3; z = Re*z;
                    end
                    if logaxis
                        q = log10(z); % plot on log scale with base
                    else
                        q = z;
                    end
                    figure;
                    if proj
                        if mesh
                            m_trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                        else
                            m_trisurf(obj.t,obj.p(:,1),obj.p(:,2),q);
                        end
                    else
                        if mesh
                            trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                        else
                            trisurf(obj.t,obj.p(:,1),obj.p(:,2),q);
                            shading flat
                        end
                        view(2);
                    end
                    cmocean('thermal',numticks(1)); cb = colorbar;
                    if logaxis
                        if length(numticks) == 3
                            desiredTicks = round(10.^(linspace(...
                                           log10(numticks(2)),...
                                           log10(numticks(3)),...
                                           numticks(1))),-1);
                        else
                            desiredTicks = round(10.^(linspace(min(q),...
                                                 max(q),numticks(1))),-1);
                        end
                        caxis([log10(min(desiredTicks)) log10(max(desiredTicks))]);
                        cb.Ticks     = log10(desiredTicks);
                        for i = 1 : length(desiredTicks)
                            cb.TickLabels{i} = num2str(desiredTicks(i));
                        end
                    elseif length(numticks) == 3
                        caxis([numticks(2) numticks(3)]);
                    end
                    ylabel(cb,'element circumradius [m]','fontsize',15);
                    title('mesh resolution');
                case('resodx')
                    TR = triangulation(obj.t,obj.p(:,1),obj.p(:,2));
                    [cc,cr] = circumcenter(TR);
                    
                    for i = 1  : length(obj.t)
                        cl = cr(i) ;
                        for j = 1 :  3
                            z(obj.t(i,j),1) = cl;
                        end
                    end
                    [Hx,Hy] = Unstruc_Bath_Slope( obj.t,obj.p(:,1),obj.p(:,2),z);
                    HH = sqrt(Hx.^2 + Hy.^2);
                    if proj
                        figure; m_trimesh(obj.t,obj.p(:,1),obj.p(:,2),HH);
                    else
                        figure; trimesh(obj.t,obj.p(:,1),obj.p(:,2),HH,...
                            'facecolor', 'flat', 'edgecolor', 'none');
                    end
                    cmocean('balance',5);
                    caxis([0 0.25]);
                    cb = colorbar;
                    title('Relaxation rate of topology');
                    ylabel(cb,'decimal percent');
                    
                case('tau0')
                    if ~isempty(obj.f13)
                        ii = find(contains({obj.f13.defval.Atr(:).AttrName},'primitive'));
                        defval  = obj.f13.defval.Atr(ii).Val;
                        userval = obj.f13.userval.Atr(ii).Val;
                        figure;
                        fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),userval(2,:)');
                        colormap([1 0 0; 0 1 0; 0 0 1]);
                        colorbar;
                    else
                        display('Fort13 structure is empty!');
                    end
                case('itfric')
                    if ~isempty(obj.f13)
                        ii = find(contains({obj.f13.defval.Atr(:).AttrName},'internal'));
                        userval = obj.f13.userval.Atr(ii).Val;
                        values = max(userval(2:end,:)',[],2);
                        [~,bpt] = extdom_edges2(obj.t,obj.p);
                        figure;
                        if proj
                            m_fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),values);
                            hold on, m_plot(bpt(:,1),bpt(:,2),'r.');
                        else
                            fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),values);
                            hold on, plot(bpt(:,1),bpt(:,2),'r.');
                        end
                        colormap(cmocean('deep',10));
                        caxis([0 5e-5])
                        colorbar;
                    end
                case('cfvals')
                    if ~isempty(obj.f13)
                        ii = find(contains({obj.f13.defval.Atr(:).AttrName},'quadratic'));
                        defval  = obj.f13.defval.Atr(ii).Val;
                        userval = obj.f13.userval.Atr(ii).Val;
                        values = obj.p(:,1)*0 + defval;
                        values(userval(1,:)) = userval(2,:);
                        figure;
                        fastscatter(obj.p(:,1),obj.p(:,2),values);
                        nouq = length(unique(values));
                        colormap(jet(nouq));
                        colorbar;
                    else
                        display('Fort13 structure is empty!');
                    end
                case('mann')
                    if ~isempty(obj.f13)
                        ii = find(contains({obj.f13.defval.Atr(:).AttrName},'mann'));
                        defval  = obj.f13.defval.Atr(ii).Val;
                        userval = obj.f13.userval.Atr(ii).Val;
                        values = max(userval(2:end,:)',[],2);
                        alltogether = zeros(np_g,1)+defval ;
                        alltogether(userval(1,:)',1) = values;
                        if proj
                            figure;
                            m_trisurf(obj.t,obj.p(:,1),obj.p(:,2),alltogether(kept));
                        else
                            trisurf(obj.t,obj.p(:,1),obj.p(:,2),alltogether(kept));
                        end
                        nouq = length(unique(values));
                        colormap(jet(nouq));
                        colorbar;
                        title('Manning n')
                    else
                        display('Fort13 structure is empty!');
                    end
                case('sponge')
                    ii = find(contains({obj.f13.defval.Atr(:).AttrName},'sponge'));
                    defval  = obj.f13.defval.Atr(ii).Val;
                    userval = obj.f13.userval.Atr(ii).Val;
                    values = userval(2,:);
                    figure;
                    fastscatter(obj.p(:,1),obj.p(:,2),defval(1)*ones(length(obj.p),1));
                    hold on
                    fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),values');
                    colormap(cmocean('deep'));
                    %                     caxis([0 5e-5])
                    colorbar;
                case('transect')
                    if proj
                        error('To plot transects, you must plot with proj=0!');
                    end
                    cmap = cmocean('ice',256);
                    figure;
                    subplot(2,1,1)
                    trisurf(obj.t,obj.p(:,1),obj.p(:,2),obj.b); view(2);
                    alpha 0.75
                    shading interp;
                    colormap(cmap); cb=colorbar; ylabel(cb,'m below geoid');
                    axis equal
                    h = imline;
                    pos = h.getPosition;
                    hold on; text(pos(1,1),pos(1,2),'Start','color','r');
                    hold on; text(pos(2,1),pos(2,2),'End','color','r');
                    [la,lo]=interpm(pos(:,2),pos(:,1),5/111e3);
                    hold on; plot(lo,la,'r.')
                    transect = [lo,la]; clearvars lo la
                    if ~isempty(obj.b)
                        F = scatteredInterpolant(obj.p(:,1),obj.p(:,2),obj.b);
                    else
                        error('Bathy must be on mesh to plot transects!');
                    end
                    bOnTransect = F(transect) ;
                    subplot(2,1,2);
                    plot(-bOnTransect,'linewi',2) ;
                    disp(mean(bOnTransect))
                    title('Bathymetry along transect');
                    ylabel('m below geoid');
                    xlabel('Points along transect')
                otherwise
                    error('Specified type is incorrect');
            end
            if proj == 1
                % now add the box
                m_grid('FontSize',16); %'box','none') %,'FontSize',12);
            end
        end
        
        % renumber mesh min. bw
        function obj=renum(obj)
            
            np=length(obj.p); nt=length(obj.t);
            % Calculate adjacency matrix of t
            S = sparse(obj.t(:,[1,1,2,2,3,3]),obj.t(:,[2,3,1,3,1,2]),1,np,np);
            W = sum(S,2);
            if any(W==0)
                error('Invalid mesh. Hanging nodes found. Retriangulate.');
            end
            % calc bw
            [i,j] = find(S);
            bw = max(i-j) + 1;
            disp(['Initial bandwidth is ',num2str(bw)]);
            
            % renumber with minimizing bw
            perm = symrcm(S);
            R    = S(perm,perm);
            prn  = obj.p(perm,:);
            if ~isempty(obj.b)
                brn  = obj.b(perm,:);
                obj.b    = brn;
            end
            
            if ~isempty(obj.bx)
                brn  = obj.bx(perm,:);
                obj.bx    = brn;
                
                brn  = obj.by(perm,:);
                obj.by    = brn;
            end
            obj.p    = prn;
            perm_inv(perm(1:np)) = ( 1 : np );
            ttemp = obj.t ;
            for ie = 1 : nt
                nm1 = ttemp(ie,1) ;
                nm2 = ttemp(ie,2) ;
                nm3 = ttemp(ie,3) ;
                ttemp(ie,1) = perm_inv(nm1);
                ttemp(ie,2) = perm_inv(nm2);
                ttemp(ie,3) = perm_inv(nm3);
            end
            obj.t = ttemp ;
            % compare bw
            [i,j] = find(R);
            bw = max(i-j) + 1;
            disp(['Renumbered bandwidth is ',num2str(bw)]);
            
            % renum the nodestrings
            if ~isempty(obj.op)
                disp('Renumbering the elevation specified boundary nodestrings...');
                for ib = 1 : obj.op.nope
                    %for iv = 1 : obj.op.nvdll(ib)
                    node = obj.op.nbdv(1:obj.op.nvdll(ib),ib);
                    obj.op.nbdv(1:obj.op.nvdll(ib),ib) = perm_inv(node);
                    %end
                end
            end
            
            if ~isempty(obj.bd)
                disp('Renumbering the flux boundary nodestrings...');
                for ib = 1 : obj.bd.nbou
                    tempcell{ib} = obj.bd.nbvv(1:obj.bd.nvell(ib),ib);
                end
                temp = cell2mat(tempcell');
                ex   = find(temp~=0) ;
                temp(ex) = perm_inv(temp(ex))';
                temp = mat2cell(temp,cellfun(@length,tempcell));
                nbvv = zeros(size(obj.bd.nbvv,1),size(obj.bd.nbvv,2));
                for ib = 1 : obj.bd.nbou
                    for iv = 1 : obj.bd.nvell(ib)
                        nbvv(iv,ib) = temp{ib}(iv,:);
                    end
                end
                obj.bd.nbvv = nbvv;
                
                % renumber the various components of the weir connectivity
                if isfield(obj.bd,'ibconn')
                    for ib = 1 : obj.bd.nbou
                        tempcell{ib} = obj.bd.ibconn(1:obj.bd.nvell(ib),ib);
                    end
                    temp = cell2mat(tempcell');
                    ex   = find(temp~=0) ;
                    temp(ex) = perm_inv(temp(ex))';
                    temp = mat2cell(temp,cellfun(@length,tempcell));
                    ibconn = zeros(size(obj.bd.nbvv,1),size(obj.bd.nbvv,2));
                    for ib = 1 : obj.bd.nbou
                        for iv = 1 : obj.bd.nvell(ib)
                            ibconn(iv,ib) = temp{ib}(iv,:);
                        end
                    end
                    obj.bd.ibconn = ibconn;
                end
                
                if isfield(obj.bd,'barinht')
                    for ib = 1 : obj.bd.nbou
                        tempcell{ib} = obj.bd.barinht(1:obj.bd.nvell(ib),ib);
                    end
                    temp = cell2mat(tempcell');
                    ex   = find(temp~=0) ;
                    temp(ex) = perm_inv(temp(ex))';
                    temp = mat2cell(temp,cellfun(@length,tempcell));
                    barinht = zeros(size(obj.bd.nbvv,1),size(obj.bd.nbvv,2));
                    for ib = 1 : obj.bd.nbou
                        for iv = 1 : obj.bd.nvell(ib)
                            barinht(iv,ib) = temp{ib}(iv,:);
                        end
                    end
                    obj.bd.barinht = barinht;
                end
                
                if isfield(obj.bd,'barincfsb')
                    for ib = 1 : obj.bd.nbou
                        tempcell{ib} = obj.bd.barincfsb(1:obj.bd.nvell(ib),ib);
                    end
                    temp = cell2mat(tempcell');
                    ex   = find(temp~=0) ;
                    temp(ex) = perm_inv(temp(ex))';
                    temp = mat2cell(temp,cellfun(@length,tempcell));
                    barincfsb = zeros(size(obj.bd.nbvv,1),size(obj.bd.nbvv,2));
                    for ib = 1 : obj.bd.nbou
                        for iv = 1 : obj.bd.nvell(ib)
                            barincfsb(iv,ib) = temp{ib}(iv,:);
                        end
                    end
                    obj.bd.barincfsb = barincfsb;
                end
                
                if isfield(obj.bd,'barincfsp')
                    for ib = 1 : obj.bd.nbou
                        tempcell{ib} = obj.bd.barincfsp(1:obj.bd.nvell(ib),ib);
                    end
                    temp = cell2mat(tempcell');
                    ex   = find(temp~=0) ;
                    temp(ex) = perm_inv(temp(ex))';
                    temp = mat2cell(temp,cellfun(@length,tempcell));
                    barincfsp = zeros(size(obj.bd.nbvv,1),size(obj.bd.nbvv,2));
                    for ib = 1 : obj.bd.nbou
                        for iv = 1 : obj.bd.nvell(ib)
                            barincfsp(iv,ib) = temp{ib}(iv,:);
                        end
                    end
                    obj.bd.barincfsb = barincfsp;
                end
                
            end
            
            if ~isempty(obj.f13)
                disp('Renumbering the fort.13...');
                for i = 1 : obj.f13.nAttr
                    idx=obj.f13.userval.Atr(i).Val(1,:);
                    idx=perm_inv(idx);
                    obj.f13.userval.Atr(i).Val(1,:) = idx;
                end
            end
        end
        
        % interp bathy/slope
        function obj = interp(obj,geodata,varargin)
            % if give cell of geodata or dems then interpolate all
            if iscell(geodata) || isstring(geodata)
                for i = 1:length(geodata)
                    if isempty(varargin)
                        obj = GridData(geodata{i},obj);
                    else
                        obj = GridData(geodata{i},obj,varargin);
                    end
                end
            else
                if isempty(varargin)
                    obj = GridData(geodata,obj);
                else
                    obj = GridData(geodata,obj,varargin);
                end
            end
        end
        
        function [obj,qual] = clean(obj,varargin)
            % [obj,qual] = clean(obj,varargin)
            % obj - msh object
            % varargin - optional base cleaning type followed by optional
            % name-value pairs as listed below:
            %
            % base cleaning types: 'medium' (or 'default'), 'passive', 'aggressive'
            %
            % optional name-value pairs
            % 'db'  - boundary element cutoff quality (0 - 1)
            % 'ds'  - perform direct smoother? (0 or 1)
            % 'con' - upper bound on connectivity (6-19)
            % 'djc' - dj_cutoff (0 - 1 [area portion] or > 1 [km^2])
            % 'sc_maxit' - max iterations for deletion of singly connected
            %         elements ( >= 0, if set to 0 operation not performed)
            % 'mqa' - allowable minimum element quality (0 - 1); setting
            %         this too value high may prevent convergence
            % 'nscreen' - print info to screen? (default = 1)
            % 'pfix' - fixed points to keep (default empty)
            % 'proj' -to project or not (default = 1)
            
            if ~isempty(varargin)
                if iscell(varargin{1}); varargin = varargin{1}; end
            end
            % keep for parsing to recursive function
            varargino = varargin;
            
            % Fixing up the mesh automatically
            disp('Beginning mesh cleaning and smoothing operations...');
 
            %process categorical cleaning options
            if any(strcmp(varargin,'passive'))
                disp('Employing passive option')
                opt.db = 0.025; opt.ds = 0; opt.con = 10; opt.djc = 1e-4; 
                opt.sc_maxit = 0; opt.mqa = 1e-4;
                varargin(strcmp(varargin,'passive')) = [];
            elseif any(strcmp(varargin,'aggressive'))
                disp('Employing aggressive option')
                opt.db = 0.25; opt.ds = 1; opt.con = 9; opt.djc = 0.25; 
                opt.sc_maxit = inf; opt.mqa = 0.1;
                varargin(strcmp(varargin,'aggressive')) = [];
            else
                disp('Employing default (medium) option')
                opt.db = 0.1; opt.ds = 1; opt.con = 9; opt.djc = 0.1; 
                opt.sc_maxit = 1; opt.mqa = 0.025;
                varargin(strcmp(varargin,'default')) = []; 
                varargin(strcmp(varargin,'medium')) = []; 
            end
            % set defaults
            opt.nscreen = 1; opt.projL = 1; pfixV = []; 
            % process user-defined individual cleaning options
            optstring = {'nscreen','pfix','proj','con','djc','db','ds',...
                         'sc_maxit','mqa'};
            for ii = 1:2:length(varargin)
                jj = find(strcmp(varargin{ii},optstring));
                if isempty(jj)
                    error(['Unrecognized input option: ' varargin{ii}])
                elseif jj == 1
                    opt.nscreen = varargin{ii + 1};
                elseif jj == 2
                    pfixV = varargin{ii + 1};
                elseif jj == 3
                    opt.projL = varargin{ii + 1};
                elseif jj == 4
                    opt.con = varargin{ii + 1};
                elseif jj == 5
                    opt.djc  = varargin{ii + 1};
                elseif jj == 6
                    opt.db = varargin{ii + 1};
                elseif jj == 7
                    opt.ds  = varargin{ii + 1};
                elseif jj == 8
                    opt.sc_maxit  = varargin{ii + 1};
                elseif jj == 9
                    opt.mqa  = varargin{ii + 1};
                end
            end
            
            % display options
            disp('the following cleaning options have been enabled..')
            disp(opt)
            disp(['length of pfix = ' length(pfixV)])
            
            if opt.projL
                global MAP_PROJECTION MAP_VAR_LIST MAP_COORDS
                if ~isempty(obj.coord)
                    MAP_PROJECTION = obj.proj ;
                    MAP_VAR_LIST   = obj.mapvar ;
                    MAP_COORDS     = obj.coord ;
                end     
                % transform to projected coordinates 
                [obj.p(:,1),obj.p(:,2)] =  m_ll2xy(obj.p(:,1),obj.p(:,2));
                if ~isempty(pfixV)
                   [pfixV(:,1),pfixV(:,2)] = m_ll2xy(pfixV(:,1),pfixV(:,2)); 
                end
            end
            
            % "fix" mesh
            [obj.p,obj.t] = fixmesh(obj.p,obj.t);
            
            if opt.db
                LT = size(obj.t,1);
                while 1
                    % Begin by just deleting poor mesh boundary elements
                    tq = gettrimeshquan(obj.p,obj.t);
                    % Get the elements that have a boundary bar
                    bdbars = extdom_edges2(obj.t,obj.p);
                    bdnodes = unique(bdbars(:));
                    vtoe = VertToEle(obj.t);
                    bele = unique(vtoe(:,bdnodes)); bele(bele == 0) = [];
                    tqbou = tq.qm(bele);
                    % Delete those boundary elements with quality < opt.db
                    if min(tqbou) >= opt.db; break; end
                    obj.t(bele(tqbou < opt.db),:) = [];
                    [obj.p,obj.t] = fixmesh(obj.p,obj.t);
                end
                if opt.nscreen
                    disp(['Deleted ' num2str(LT-size(obj.t,1)) ...
                          ' bad boundary elements'])
                end
            end
            
            % Make mesh traversable
            obj = Make_Mesh_Boundaries_Traversable(obj,opt.djc,opt.nscreen);
            
            % Delete elements with single edge connectivity
            obj = Fix_single_connec_edge_elements(obj,opt.sc_maxit,opt.nscreen);
            
            % Reduce the mesh connectivity to maximum of con-1
            obj = renum(obj);
            % May not always work without error
            if opt.con > 6
                try
                   obj = bound_con_int(obj,opt.con);
                catch
                   warning('Could not reduce connectivity mesh');
                end
            end
           
            % Now do the smoothing if required
            if opt.ds
                % Perform the direct smoothing
                [obj.p,obj.t] = direct_smoother_lur(obj.p,obj.t,...
                                                    pfixV,opt.nscreen);
            end
           
            % Checking and displaying element quality
            tq = gettrimeshquan( obj.p, obj.t);            
            mq_m = mean(tq.qm);
            mq_l = min(tq.qm);
            mq_s = std(tq.qm);
            mq_l3sig = mq_m - 3*mq_s;
            qual = [mq_m,mq_l3sig,mq_l];
            
            if mq_l < opt.mqa
                % Need to clean it again
                disp('Poor or overlapping elements, cleaning again')
                % repeat without projecting (already projected)
                ii = find(strcmp(varargino,'proj'));
                if ~isempty(ii)
                    varargino{ii+1} = 0;
                else
                    varargino{end+1} = 'proj';
                    varargino{end+1} = 0;
                end
                obj = clean(obj,varargino(:));
            end
            
            if opt.nscreen
                disp(['number of nodes is ' num2str(length(obj.p))])
                disp(['mean quality is ' num2str(mq_m)])
                disp(['min quality is ' num2str(mq_l)])
            end
           
            % Do the transformation back
            if opt.projL
                [obj.p(:,1),obj.p(:,2)] = m_xy2ll(obj.p(:,1),obj.p(:,2));
            end
        end
        
        function obj = lim_bathy_slope(obj,dfdx,overland)
            %obj = lim_bathy_slope(obj,dfdx,overland)
            if nargin < 3
                overland  = 0;
            end
            % Limit to topo or bathymetric slope to dfdx on the edges
            imax = 100;
            [edge,elen] = GetBarLengths(obj,0);  
            bt = obj.b; 
            if overland
                I = bt > 0; 
                word = 'topographic';
            else
                I = bt < 0; 
                word = 'bathymetric';
            end
            bt(I) = 0;
            [bnew,flag] = limgrad(edge,elen,bt,dfdx,imax);
            if flag
               obj.b(~I) = bnew(~I);
               disp(['Successfully limited ' word ' slope to ' ...
                      num2str(dfdx) ' in limgrad function']) 
            else
               warning(['Could not limit ' word ' slope to ' ...
                        num2str(dfdx) ' in limgrad function']) 
            end
        end
        
                % make nodestrings
        function obj = makens(obj,type,dir,cutlim,depthlim)
            if nargin < 2
               error('Needs type: one of auto, islands, delete, or outer')
            end
            if nargin < 4 || isempty(cutlim)
              cutlim = 10 ; 
            end
            if nargin < 5 || isempty(depthlim)
              depthlim = 10 ; 
            end
            L = 1e3;
%             trim = 0; periodic = 0;
%             if strcmp(type(max(1,end-3):end),'trim')
%                 type = type(1:end-4); trim = 1;
%             end
            switch type
                case('auto')
                    if nargin < 3 || ~isa(dir,'geodata')
                        error('third input must be a geodata class for auto makens')
                    else
                        gdat = dir;
                    end
                    
                    % Check for global mesh
                    Le = find(obj.p(:,1) < -179, 1);
                    Ri = find(obj.p(:,1) > 179, 1);
                    if ~isempty(Le) && ~isempty(Ri)
                        error('Detected global mesh, cannot apply makens')
                    end
                    
                    % Get the boundaries
                    [etbv,~]  = extdom_edges2(obj.t,obj.p);
                    [~,poly_idx] = extdom_polygon(etbv,obj.p,1);
                    
                    outerbox = gdat.outer(1:find(isnan(gdat.outer(:,1)),1,'first'),:);
                    if ~isempty(gdat.mainland)
                        mainland = gdat.mainland;
                        mainland(isnan(mainland(:,1)),:) = [];
                    end
                    if ~isempty(gdat.inner)
                        inner = gdat.inner;
                        inner(isnan(inner(:,1)),:) = [];
                    end
                    
                    % Get geodata outer and mainland polygons
                    nope = 0; neta = 0; nbou  = 0; nvel  = 0;
        
                    % loop through all the polygons
                    for poly_count = 1 : length(poly_idx)
                        idv = poly_idx{poly_count};
                        % 
                        [~,odst] = ourKNNsearch(outerbox(1:end-1,:)',obj.p(idv,:)',1);         
                        if ~isempty(gdat.mainland)
                            [~,mdst] = ourKNNsearch(mainland',obj.p(idv,:)',1);
                        else
                            mdst = 1e4;
                        end                        
                        
                        if ~isempty(gdat.inner)
                            [~,idst] = ourKNNsearch(inner',obj.p(idv,:)',1);
                            if mean(idst) < mean(odst) && ...
                               mean(idst) < mean(mdst)
                                % must be an island
                                nbou = nbou + 1;
                                nvell(nbou) = length(idv);
                                nvel = nvel + nvell(nbou);
                                nbvv(1:nvell(nbou),nbou) = idv';
                                ibtype(nbou) = 21;
                                continue;
                            end
                        end 
                        if ~isempty(gdat.mainland)
                            
                            % set the bathy values
                            if ~isempty(gdat.Fb)
                                zvalue = gdat.Fb(obj.p(idv,:));
                            else
                                % Dummy all underwater
                                zvalue = 0*obj.p(idv,1) - 1e4;
                            end
                            
                            % find the boundaries that are mainland
                            main = mdst < odst | zvalue > -depthlim;
                               
                            % indices of switch
                            Cuts  = find(diff(main) ~= 0);
                    
                            % Do not include open boundary that is 
                            % smaller than cutlim vertices across
                            rm = false(length(Cuts),1);
                            if main(1); st = 1; else; st = 2; end
                            for ii = st:2:length(Cuts)
                                if ii == length(Cuts)
                                    if length(main) - Cuts(ii) + ...
                                       Cuts(1) - 1  < cutlim
                                       rm([1 end]) = 1;
                                    end
                                else
                                    if Cuts(ii+1) - Cuts(ii) < cutlim
                                        rm(ii:ii+1) = 1;
                                    end
                                end
                            end
                            Cuts(rm) = [];                              

                            if isempty(Cuts)
                                if sum(main)/length(main) > 0.5
                                    % get mainland
                                    nbou = nbou + 1;
                                    nvell(nbou) = length(idv);
                                    nvel = nvel + nvell(nbou);
                                    ibtype(nbou) = 20;
                                    nbvv(1:nvell(nbou),nbou) = idv';
                                else
                                    % get ocean
                                    nope = nope + 1;
                                    nvdll(nope) = length(idv);
                                    neta = neta + nvdll(nope);
                                    ibtypee(nope) = 0;
                                    nbdv(1:nvdll(nope),nope) = idv';
                                end
                            else
                                % loop through each cut
                                for loop = 1:length(Cuts)
                                    if loop < length(Cuts)
                                        idv_temp = idv(Cuts(loop):Cuts(loop+1));
                                    else
                                        % need to loop back to start
                                        idv_temp = [idv(Cuts(loop):end); idv(1:Cuts(1))];
                                    end
                                    % Break up idv_temp into L chunks
                                    N = ceil(length(idv_temp)/L);
                                    LE = ceil(length(idv_temp)/N);
                                    ns = 1;
                                    for nn = 1:N
                                        ne = min(ns + LE - 1,length(idv_temp));
                                        idv_t = idv_temp(ns:ne);
                                        if ~main(Cuts(loop)+1)
                                            % Get ocean
                                            nope = nope + 1;
                                            nvdll(nope) = length(idv_t);
                                            neta = neta + nvdll(nope);
                                            ibtypee(nope) = 0;
                                            nbdv(1:nvdll(nope),nope) = idv_t;
                                        else
                                            % Get mainland
                                            nbou = nbou + 1;
                                            nvell(nbou) = length(idv_t);
                                            nvel = nvel + nvell(nbou);
                                            ibtype(nbou) = 20;
                                            nbvv(1:nvell(nbou),nbou) = idv_t;
                                        end
                                        ns = ne + 1;
                                    end
                                end
                            end                                               
                        else
                            % must be ocean
                            nope = nope + 1;
                            nvdll(nope) = length(idv);
                            neta = neta + nvdll(nope);
                            ibtypee(nope) = 0;
                            nbdv(1:nvdll(nope),nope) = idv';
                        end
                    end
                    
                    if nope > 0
                        % ocean boundary
                        obj.op.nope = nope ;
                        obj.op.neta = neta ;
                        obj.op.nvdll = nvdll ;
                        obj.op.ibtype = ibtypee ;
                        obj.op.nbdv = nbdv;
                    end
                    
                    if nbou > 0
                        % land boundary
                        obj.bd.nbou = nbou ;
                        obj.bd.nvel = nvel ;
                        obj.bd.nvell = nvell ;
                        obj.bd.ibtype = ibtype ;
                        obj.bd.nbvv = nbvv ;    
                    end
                    
                case('islands')
                    [etbv,~]  = extdom_edges2(obj.t,obj.p);
                    [poly,poly_idx,max_ind] = extdom_polygon(etbv,obj.p,1);
                    
                    nbou = 0;
                    nvel = 0;
                    % the largest polygon will be a combination of ocean and mainland
                    % boundaries. Deal with this first, then remove it from the polygon
                    %poly(max_ind) = [];
                    %poly_idx(max_ind)=[];
                    
                    % loop through the remaining polygons
                    for poly_count = 1 : length(poly)
                        vso = poly{poly_count};
                        idv = poly_idx{poly_count};
                        % islands
                        nbou = nbou + 1;
                        nvell(nbou) = length(vso);
                        nvel = nvel + nvell(nbou);
                        nbvv(1:nvell(nbou),nbou) = idv';
                        ibtype(nbou) = 21;
                    end
                    
                    % ocean boundary
                    obj.op.nope = 0 ;
                    obj.op.neta = 0 ;
                    obj.op.nvdll = 0 ;
                    obj.op.ibtype = 0 ;
                    obj.op.nbdv = 0;
                    
                    % land boundary
                    if nbou == 0
                        disp('No islands found!')
                    else
                        obj.bd.nbou = nbou ;
                        obj.bd.nvel = nvel ;
                        obj.bd.nvell = nvell ;
                        obj.bd.ibtype = ibtype ;
                        obj.bd.nbvv = nbvv ;
                    end
                        
                case('delete')
                    % have the user select the nodestring '
                    plot(obj,'bd') ;
                    temp = obj.bd.nbvv;
                    bounodes=obj.bd.nbvv ;
                    idx=sum(bounodes~=0);
                    bounodes=bounodes(:) ;
                    bounodes(bounodes==0)=[] ;
                    boupts = obj.p(bounodes,:) ;
                    idx2=[0,cumsum(idx)]'+1;
                    k = 0 ;
                    for i = 1 : length(idx2)-1
                        k = k + 1 ;
                        bounodes(idx2(i):idx2(i+1)-1,2) = k ;
                    end
                    
                    dcm_obj = datacursormode(gcf);
                    title('use data cursor to select nodestring to be deleted');
                    pause
                    c_info = getCursorInfo(dcm_obj);
                    [tmp(:,1),tmp(:,2)]=m_ll2xy(boupts(:,1),boupts(:,2));
                    idx3 = ourKNNsearch(boupts',c_info.Position',1)  ;
                    del = bounodes(idx3,2) ;  %<- get the nodestring index
                    pltid = temp(:,del) ; pltid(pltid==0)=[] ;
                    hold on; m_plot(obj.p(pltid,1),obj.p(pltid,2),'r-','linewi',2) ;
                    disp(['Delete boundary with index ',num2str(del),'?']) ;
                    pause
                    obj.bd.nbvv(:,del)=[];
                    num_delnodes = idx(del) ;
                    obj.bd.nbou = obj.bd.nbou - 1 ;
                    obj.bd.nvell(del)=[] ;
                    obj.bd.ibtype(del)=[] ;
                    obj.bd.nvel = obj.bd.nvel - num_delnodes ;
                    
                case('outer')
                    if nargin < 3
                       error('must specify direction of boundary in third entry') 
                    end
                        
                    [bnde,bpts]=extdom_edges2(obj.t,obj.p);
                    
                    % use this to figure out the vstart and vend
                    figure, plot(bpts(:,1),bpts(:,2),'k.');
                    %hold on; fastscatter(obj.p(:,1),obj.p(:,2),obj.b) ;
                    caxis([-10 10]) ; axis equal ;
                    title('use data cursor to identify vstart and vend');
                    dcm_obj = datacursormode(gcf);
                    set(dcm_obj,'UpdateFcn',{@myupdatefcn2,bpts})
                    
                    % Use data cursor to get the values of the boundary nodes.
                    vstart=input('Enter the value for vstart : ');
                    vend=input('Enter the value for vend : ');
                    bndidx=unique(bnde(:));
                    
                    vstart= bndidx(vstart);
                    vend  = bndidx(vend);
                    
                    [~,~,obj.op,obj.bd] = extract_boundary(vstart,vend,bnde,obj.p,...
                        dir,obj.op,obj.bd); %<--updates op and bd.
                case('periodic')
                    % get elements that straddle the -180/180 boundary
                    bars = [obj.t(:,[1,2]); obj.t(:,[1,3]); obj.t(:,[2,3])]; % Interior bars duplicated
                    bars = unique(sort(bars,2),'rows');                  % Bars as node pairs
                    dlon = obj.p(bars(:,1),1)- obj.p(bars(:,2),1);     % List of bar vectors
                    % Only add the periodic point list
                    %[I,d] = ourKNNsearch(obj.p',obj.p'+[360;0],1);
                    %J = find(d < 1e-5);
                    % list of points that straddle the -180/180 boundary
                    periodic_bc_list = bars(abs(dlon) > 180,:);
                    obj.bd.nbou = 1 ;
                    obj.bd.nvel = length(periodic_bc_list) ;
                    obj.bd.nvell = obj.bd.nvel ;
                    obj.bd.ibtype = 94 ;
                    obj.bd.nbvv = periodic_bc_list ;
                    return;
                    
                case('weirs')
                    if ~isa(dir,'geodata')
                        error('The third input must be a geodata class object you used create the mesh with.')
                    else
                        gdat = dir;
                    end
                    
                    % identifying and adding internal weir type boundaries (ibtype=24)
                    for ii = 1 : length(gdat.ibconn_pts) % for each weir
                        [front_nn, d1] = ourKNNsearch(obj.p',gdat.ibconn_pts{ii}(:,1:2)',1);
                        [back_nn,  d2] = ourKNNsearch(obj.p',gdat.ibconn_pts{ii}(:,3:4)',1);
                        rm = d1 > 1e-9 & d2 > 1e-9;
                        front_nn(rm) = [] ; back_nn(rm) = [] ;
                        nn = [front_nn,back_nn] ;
                        for iii = 1 : length(nn)
                            rm2(iii,1)=length(unique(nn(iii,:)))~=2 ;
                        end
                        front_nn(rm2) = [] ; back_nn(rm2) = [] ;
                        clearvars rm2 ;
                        %                         % visualize node pairs
                        %                         plot(obj,'tri',0) ;
                        %                         hold on; plot(obj.p(front_nn,1),obj.p(front_nn,2),'r.') ;
                        %                         hold on; plot(obj.p(back_nn,1),obj.p(back_nn,2),'gs') ;
                        
                        % add to the boundary struct
                        if  ~isempty(obj.bd)
                            % unpack current boundaries
                            nbou   = obj.bd.nbou;
                            nvel   = obj.bd.nvel;
                            nvell  = obj.bd.nvell;
                            nbvv   = obj.bd.nbvv;
                            ibtype = obj.bd.ibtype;
                            
                            % append on the weir
                            nbou = nbou + 1;
                            nvell(nbou) = length(front_nn) ; % only list front facing nodes in nbvv for barriers
                            nvel = nvel + nvell(nbou);
                            nbvv(1:nvell(nbou),nbou) = front_nn;
                            ibtype(nbou) = 24;
                            
                            % repackage it
                            obj.bd.nbou   = nbou ;
                            obj.bd.nvel   = nvel ;
                            obj.bd.nvell  = nvell ;
                            obj.bd.ibtype = ibtype ;
                            obj.bd.nbvv   = nbvv ;
                        else
                            % create a new boundary struct
                            % initialize arrays
                            nbou  = 0 ;
                            nvel  = 0 ;
                            nvell = [] ;
                            ibtype= [] ;
                            nbvv  = [] ;
                            
                            % append on the weir
                            nbou = nbou + 1;
                            nvell(nbou) = length(front_nn) ; % only list front facing nodes in nbvv for barriers
                            nvel = nvel + nvell(nbou);
                            nbvv(1:nvell(nbou),nbou) = front_nn;
                            ibtype(nbou) = 24;
                            
                            % package it
                            obj.bd.nbou   = nbou ;
                            obj.bd.nvel   = nvel ;
                            obj.bd.nvell  = nvell ;
                            obj.bd.ibtype = ibtype ;
                            obj.bd.nbvv   = nbvv ;
                        end
                        
                        % either append or create new weir connectivity tables
                        if isfield(obj.bd,'ibconn')
                            % then append
                            ibconn    = obj.bd.ibconn ;
                            barinht   = obj.bd.barinht ;
                            barincfsb = obj.bd.barincfsb ;
                            barincfsp = obj.bd.barincfsp ;
                            ibconn(1:nvell(nbou),nbou) = back_nn ; % only list back facing nodes in ibconn for barriers
                            % Opt 1) ask the user for a dataset to give the
                            % crestline (barinht) ;
                            ar = input('Type 1 to enter weir crest height or type 2 to specify dataset...') ;
                            if ar==1
                                % fixed value for weir crests
                                ht = input('Enter value in meters ABOVE the geoid for the height of the weir...') ;
                                disp('-----------------------------------------------------------') ;
                            elseif ar==2
                                % working on it !
                                %
                                error('NOT WORKING YET')
                            end
                            barinht(1:nvell(nbou),nbou) =  ht ;
                            barincfsb(1:nvell(nbou),nbou) = 1 ; % these are standard values
                            barincfsp(1:nvell(nbou),nbou) = 1 ; % these are standard values
                        else
                            % then create new
                            ibconn = nbvv*0 ;
                            ibconn(1:nvell(nbou),nbou) = back_nn ;  % only list back facing nodes in ibconn for barriers
                            barinht =  [] ;
                            barincfsb = [];
                            barincfsp = [] ;
                            % Opt 1) ask the user for a dataset to give the crestline
                            % height (barinht)
                            ar = input('Type 1 to enter weir crest height or type 2 to specify dataset...') ;
                            if ar==1
                                % fixed value for weir crests
                                ht = input('Enter value in meters ABOVE the geoid for the height of the weir...') ;
                                disp('-----------------------------------------------------------') ;
                            elseif ar==2
                                % working on it !
                                %
                                error('NOT WORKING YET')
                            end
                            barinht(1:nvell(nbou),nbou)   = ht ;
                            barincfsb(1:nvell(nbou),nbou) = 1 ; % these are standard values
                            barincfsp(1:nvell(nbou),nbou) = 1 ; % these are standard values
                        end
                        % now add them back
                        obj.bd.ibconn    = ibconn ;
                        obj.bd.barinht   = barinht ;
                        obj.bd.barincfsb = barincfsb ;
                        obj.bd.barincfsp = barincfsp ;
                    end
            end
            
            
            function txt = myupdatefcn2(~,event_obj,myarray)
                pos = get(event_obj,'Position');
                ind = find(abs(myarray(:,1)-pos(1))<eps & abs(myarray(:,2)-pos(2))<eps);
                txt = {['X: ',num2str(pos(1))],...
                    ['Y: ',num2str(pos(2))],...
                    ['Index: ',num2str(ind')]};
            end
        end
        
   
        
        function merge = plus(obj1,obj2,tight)
            % merge = plus(obj1,obj2,tight)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Merge together two meshes contained in the msh objects obj1
            % and obj2. It uses MATLAB's implementation of the Boywer-Watson
            % incremental triangulation and then applies mesh cleaning
            % algorithms to "fix" the intersection zone between meshes.
            %
            % NOTE: msh obj1 can contain fixed points and these will be
            % preserved in the merger.
            %
            % HINT: if you intend to use this a priori, perhaps you should
            % ensure the mesh size functions are similar.
            %
            % INPUTS:
            % obj1: msh() class of inset mesh (INSET MUST BE FIRST).
            % obj2: msh() class of base mesh.
            % tight: = 1 (default) updates boundary of obj2 (for removal of 
            %          triangles) after removing intersection(obj1,obj2) 
            %          from obj2 (best for merging overlapping meshes that 
            %          have similar resolution/mesh size functions)
            %        = 0 does not update boundary of obj2 (for removal of 
            %          triangles) after removing insersection(obj1,obj2) 
            %          from obj2 (good for avoiding holes when merging 
            %          meshes with disparate resolution/mesh size functions)
            %        = -1 does not remove intersection(obj1,obj2) from obj2
            %          (assumes there is no overlap, best for merging
            %          [almost] non-overlapping meshes - ought to have 
            %          matching or close to matching vertices but this is 
            %          NOT a requirement)
            %
            % OUPUTS:
            % merge: a msh object in which the msh obj1's connectivity and 
            % bathymetry is carried over (along with fixed points and edges
            % from obj1).
            %
            % kjr, und, chl, sept. 2017 Version 1.0.
            % kjr, und, chl, oct. 2018, Version 1.5
            % kjr, usp. july 2019. Support for carrying over weirs and pfix. 
            % wjp, inserting the 'tight' options optimizied for certain
            % situations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if nargin < 3
                tight = 1;
            end
            
            p1 = obj1.p; t1 = obj1.t;
            p2 = obj2.p; t2 = obj2.t;
            
            pfix1 = obj1.pfix; nfix1 = length(pfix1); egfix1 = obj1.egfix; 
            pfix2 = obj2.pfix; nfix2 = length(pfix2); egfix2 = obj2.egfix; 
            
            % all combined edge constraints 
            nfix   = nfix1+nfix2 ; 
            egfixx  = [egfix1; egfix2+nfix1] ; 
            pfixx   = [pfix1; pfix2] ;
            
            global MAP_PROJECTION MAP_COORDS MAP_VAR_LIST
            if ~isempty(obj2.coord)
                if any(min(obj1.p) < min(obj2.p)) || ...
                        any(max(obj1.p) > max(obj2.p))
                    objt = obj2; objt.p = [objt.p; obj1.p];
                    objt.p(1,:) = min(objt.p) - 1e-3;
                    objt.p(2,:) = max(objt.p) + 1e-3;
                    setProj(objt,1,obj2.proj.name);
                else
                    % kjr 2018,10,17; Set up projected space imported from msh class
                    MAP_PROJECTION = obj2.proj ;
                    MAP_VAR_LIST   = obj2.mapvar ;
                    MAP_COORDS     = obj2.coord ;
                end
            else
                projname = 'stereo';
                setProj(obj2,1,projname);
            end
           
            % project both meshes into the space of the global mesh
            [p1(:,1),p1(:,2)] = m_ll2xy(p1(:,1),p1(:,2)) ;
            [p2(:,1),p2(:,2)] = m_ll2xy(p2(:,1),p2(:,2)) ;
            if nfix > 0
                [pfixx(:,1),pfixx(:,2)] = m_ll2xy(pfixx(:,1),pfixx(:,2)) ;
            end
            
            [~,d1] = ourKNNsearch(p1',p1',2);
            [~,d2] = ourKNNsearch(p2',p2',2);
            mvd1 = max(d1); mvd2 = max(d2);
            midi = min(min(d2),min(d1)); 
            overlap = 0.5*(mvd1 + mvd2);
            
            disp('Forming outer boundary for base...')
            try
                cell2 = extdom_polygon(extdom_edges2(t2,p2),p2,-1,0);
                poly_vec2 = cell2mat(cell2'); 
                [edges2] = Get_poly_edges(poly_vec2);
            catch
                error('Mesh 2 is invalid. Please execute msh.clean on object');
            end
            
            
            disp('Forming outer boundary for inset...')
            try
                cell1 = extdom_polygon(extdom_edges2(t1,p1),p1,-1,0);
                poly_vec1 = cell2mat(cell1');
                [edges1] = Get_poly_edges(poly_vec1);
            catch
                error('Mesh 1 is invalid. Please execute msh.clean on object');
            end
            
            % Delete the region in the global mesh that is in the
            % intersection with inset.
            if tight > -1
                disp('Calculating intersection...');
                [x3,y3] = polybool('intersection',...
                    poly_vec1(:,1),poly_vec1(:,2),poly_vec2(:,1),poly_vec2(:,2));
                if ~isempty(x3)
                    poly_vec3 = [x3,y3]; poly_vec3(end+1,:) = NaN;
                    [edges3]  = Get_poly_edges(poly_vec3);
                    in1 = inpoly(p2(t2(:,1),:),poly_vec3,edges3);
                    in2 = inpoly(p2(t2(:,2),:),poly_vec3,edges3);
                    in3 = inpoly(p2(t2(:,3),:),poly_vec3,edges3);
                    t2(in1 & in2 & in3,:) = [];
                    % We need to delete straggling elements that are
                    %  generated through the above deletion step
                    pruned2 = msh() ; pruned2.p = p2; pruned2.t = t2;
                    pruned2 = Make_Mesh_Boundaries_Traversable(pruned2,0.25,1);
                    t2 = pruned2.t; p2 = pruned2.p;
                    % get new poly_vec2
                    if tight
                        cell2 = extdom_polygon(extdom_edges2(t2,p2),p2,-1,0);
                        poly_vec2 = cell2mat(cell2');
                        [edges2] = Get_poly_edges(poly_vec2);
                    end
                end
            end
            pv1 = poly_vec1(~isnan(poly_vec1(:,1)),:);
            pv2 = poly_vec2(~isnan(poly_vec2(:,1)),:);
            
            disp('Merging...')
            if isempty(egfix1)
                DTbase = delaunayTriangulation(p1);
            else
                DTbase = delaunayTriangulation(p1,egfix1);
            end
            DTbase.Points(end+(1:length(p2)),:) = p2;

            % delete points that are too close together
            [kk,dst] = ourKNNsearch(DTbase.Points',DTbase.Points',2);
            dst = dst(:,2); kk = kk(:,2);
            del = kk(dst < midi); 
            DTbase.Points(double(unique(del)),:) = [];
            
            % ensure obj2's edge constraints are obeyed.
            if ~isempty(egfix2)
                DTbase = delaunayTriangulation(DTbase.Points,egfix2+length(p1));
                % we must shuffle around the Points in DTbase to ensure all
                % pfix are preprended at the top of the Points array.
                pm = DTbase.Points ; 
                tm = DTbase.ConnectivityList ; 
                pm = fixmesh([pfixx; pm]) ; % this will delete points that are duplicate later on in the list
                DTbase = delaunayTriangulation(pm,egfixx);
            end
            
            tq.qm = 0;
            while min(tq.qm) < 1e-4
                % Prune triangles outside both domains.
                disp('Pruning...')
                
                for ii = 1:2
                    % The loop makes sure to remove only small connectivity for the boundaries
                    if ii == 2
                        % To remove small connectivity or bad elements
                        [~, enum] = VertToEle(tm);
                        bdbars = extdom_edges2(tm,pm);
                        bdnodes = unique(bdbars(:));
                        I = find(enum <= 4);
                        nn = setdiff(I',[(1:nfix)';bdnodes]);  
                        DTbase.Points(unique(nn),:) = [];
                    end
                    
                    pm = DTbase.Points; tm = DTbase.ConnectivityList;
                    nei = DTbase.neighbors;
                    
                    pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/3;
                    
                    %in1 is inside the inset boundary polygon
                    in1 = inpoly(pmid,poly_vec1,edges1);
                    
                    %in2 is inside the global boundary polygon
                    in2 = inpoly(pmid,poly_vec2,edges2);
                    
                    %in3 is inside the intersection
                    if exist('poly_vec3','var') && tight == 1
                        in3 = inpoly(pmid,poly_vec3,edges3);
                    else
                        in3 = false(size(in1)); 
                    end
                    
                    % remove triangles that aren't in the global mesh or 
                    % aren't in the inset mesh
                    del = (~in1 & ~in2) | (~in1 & in3);
                    tm(del,:) = [];
                end
                
                merge = msh() ; merge.p = pm; merge.t = tm ;
				merge = clean(merge,'passive','proj',0,'pfix',pfixx);
                
                % if we don't know the overlap region
                [~,dst1] = ourKNNsearch(p1',merge.p',1);
                [~,dst2] = ourKNNsearch(p2',merge.p',1);
                idx = find(abs(dst1 - dst2) < overlap);
                
                % use smoother around intersection while obeying constraints)
                constr = setdiff((1:length(merge.p))',idx);
                if nfix > 0
                   constr = [(1:nfix)'; constr] ; 
                   constr = unique(constr);
                end
                [pm,tm] = smoothmesh(merge.p,merge.t,constr,50,0.01);
                tq = gettrimeshquan(pm,tm);
                disp(['min element quality is ', num2str(min(tq.qm))])
                if min(tq.qm) < 1e-4
                    DTbase = delaunayTriangulation(pm(:,1),pm(:,2));
                end
            end
            
            merge.p = pm; merge.t = tm;
            merge.pfix  = pfixx ; merge.egfix = egfixx ;
            
            % convert back to lat-lon wgs84
            [merge.p(:,1),merge.p(:,2)] = ...
                m_xy2ll(merge.p(:,1),merge.p(:,2));
            if nfix > 0
                [merge.pfix(:,1),merge.pfix(:,2)] = ...
                    m_xy2ll(pfixx(:,1),pfixx(:,2));
            end
            
            merge.proj    = MAP_PROJECTION ;
            merge.coord   = MAP_COORDS ;
            merge.mapvar  = MAP_VAR_LIST ;                        
                                    
            % Check element order
            merge = CheckElementOrder(merge);             
            
            % Carry over bathy and gradients
            merge.b = 0*merge.p(:,1);
            
            [idx1,dst1] = ourKNNsearch(obj1.p',merge.p',1);   
            [idx2,dst2] = ourKNNsearch(obj2.p',merge.p',1);   
            merge.b( dst1 <= dst2) = obj1.b( idx1(dst1 <= dst2)); 
            merge.b( dst2 <  dst1 ) = obj2.b( idx2(dst2 <  dst1) );
            % ensure depth of the first mesh is preserved as a priority
            in1 = inpoly(pm,poly_vec1,edges1);
            merge.b( in1 ) = obj1.b( idx1(in1) ); 
            if ~isempty(obj1.bx) && ~isempty(obj2.bx)
                merge.bx = 0*merge.b; merge.by = 0*merge.b;
                merge.bx(dst2 < dst1) = obj2.bx(idx2(dst2 < dst1)); 
                merge.bx(dst1 <= dst2) = obj1.bx(idx1(dst1 <= dst2)); 
                merge.bx( in1 ) = obj1.bx( idx1(in1) ); 
                merge.by(dst2 < dst1) = obj2.by(idx2(dst2 < dst1)); 
                merge.by(dst1 <= dst2) = obj1.by(idx1(dst1 <= dst2)); 
                merge.by( in1 ) = obj1.by( idx1(in1) ); 
            end
            disp(['Note that f13, f15 and boundary conditions etc. have' ...
                  'not been carried over into the merged mesh'])                                  
        end
      
        function obj = CheckElementOrder(obj,proj)
            if nargin == 1
               proj = 0; 
            end
            vx = obj.p(:,1); vy = obj.p(:,2);
            xt = [vx(obj.t(:,1)) vx(obj.t(:,2)) ...
                  vx(obj.t(:,3)) vx(obj.t(:,1))];
            yt = [vy(obj.t(:,1)) vy(obj.t(:,2)) ...
                  vy(obj.t(:,3)) vy(obj.t(:,1))];
            dxt = diff(xt,[],2);
            dyt = diff(yt,[],2);
            for ii = 1:3
                iii = ii - 1; if iii == 0; iii = 3; end
                I = abs(dxt(:,ii)) > 180 & abs(dxt(:,iii)) > 180;
                xt(I,ii) = xt(I,ii) + 360*sign(dxt(I,ii)); 
            end
            xt(:,end) = xt(:,1);
            % Get new diff
            dxt = diff(xt,[],2);
            Area = dxt(:,3).*-dyt(:,2) + dxt(:,2).*dyt(:,3);
            if ~isempty(find(Area < 0, 1))
                disp([num2str(sum(Area < 0)) ...
                 ' elements found that are not in counterclockwise order'])
                if proj
                    global MAP_PROJECTION MAP_VAR_LIST MAP_COORDS
                    % kjr 2018,10,17; Set up projected space imported from msh class
                    MAP_PROJECTION = obj.proj ;
                    MAP_VAR_LIST   = obj.mapvar ;
                    MAP_COORDS     = obj.coord ;
                    [pt(:,1),pt(:,2)] = m_ll2xy(obj.p(:,1),obj.p(:,2));
                else
                    pt = obj.p;
                end
                [~,tt] = fixmesh(pt,obj.t);
                obj.t(Area < 0,:) = tt(Area < 0,:);
                obj = CheckElementOrder(obj,~proj);
            else
                disp('All elements now in counterclockwise order')
            end
        end
            
        function [out1,barlen,bars] = CalcCFL(obj,dt,type)
            if nargin < 3
                % use spherical haversine distances
                type = 0;
            end
            g      = 9.81;        % gravity
            [bars,barlen] = GetBarLengths(obj,type);
            % sort bar lengths in ascending order
            [barlen,IA] = sort(barlen,'ascend');
            bars = bars(IA,:);
            % get the minimum bar length for each node
            [B1,IB] = unique(bars(:,1),'first');
            [B2,IC] = unique(bars(:,2),'first');
            d1 = NaN*obj.p(:,1); d2 = NaN*obj.p(:,1);
            d1(B1) = barlen(IB); d2(B2) = barlen(IC);
            d = min(d1,d2);
            
            % wavespeed in ocean (second term represents orbital
            % velocity at 0 degree phase for 1-m amp. wave).
            U = sqrt(g*max(obj.b,1)) + sqrt(g./max(obj.b,1));
            if nargin > 1 && ~isempty(dt)
                % Get CFL from input dt
                CFL = dt*U(:)./d;  % <-- from the wave celerity
                out1 = CFL;
            else
                CFL = 1.0;
                dt = d.*CFL./U; % <-- from the wave celerity
                out1 = dt;
            end
        end
        
        function [bars,barlen] = GetBarLengths(obj,type)
            % [bars,barlen] = GetBarLengths(obj,type)
            % Get bars and bar lengths of the elements in msh object
            % set type = 0 for bar lengths computed using Harvesine formula
            % set type = 1 for bar lengths computed using CCP with
            % correction factor for x-direction
            if nargin < 2
               type = 0; 
            end
            bars = [obj.t(:,[1,2]); obj.t(:,[1,3]); obj.t(:,[2,3])]; % Interior bars duplicated
            bars = unique(sort(bars,2),'rows');                      % Bars as node pairs
            if type == 0
                % Compute based on Haversine spherical earth distances
                long   = zeros(length(bars)*2,1);
                lat    = zeros(length(bars)*2,1);
                long(1:2:end) = obj.p(bars(:,1),1);
                long(2:2:end) = obj.p(bars(:,2),1);
                lat(1:2:end)  = obj.p(bars(:,1),2);
                lat(2:2:end)  = obj.p(bars(:,2),2);
                % Get spherical earth distances for bars
                barlen = m_lldist(long,lat); 
                barlen = barlen(1:2:end)*1e3;  % L = Bar lengths in meters
            elseif type == 1
                % Compute based on CPP with sfac correction factor for x-direction
                sfea = obj.p(:,2); sfac = cosd(mean(sfea)) ./ cosd(sfea);
                sfacelemid = mean(sfac(obj.t),2);
                vtoe = VertToEle(obj.t);
                vtoe(vtoe == 0) = length(obj.t) + 1;
                sfacelemid(end+1) = NaN;
                conbar1 = vtoe(:,bars(:,1));
                conbar2 = vtoe(:,bars(:,2));
                sfacbar1 = max(sfacelemid(conbar1));
                sfacbar2 = max(sfacelemid(conbar2));
                sfac = max(sfacbar1,sfacbar2)';
                % add safety factor for high latitudes (this is quite critical)
                sfac = sfac.*ceil(sfac/10); 
                [x, y] = CPP_conv( obj.p(:,1), obj.p(:,2) );
                barlen = hypot(diff(x(bars),[],2)./sfac,diff(y(bars),[],2));
            end
            
            function [ x,y ] = CPP_conv( lon,lat )
                %CPP_conv Converts lat and lon to x and y
                lon0 =  mean(lon) * pi/180; lat0 = mean(lat) * pi/180;

                R = 6378206.4;
                lonr = lon * pi/180; latr = lat * pi/180; 
                x = R * (lonr - lon0) * cos(lat0);
                y = R * latr;
            end
        end

        function obj = CheckTimestep(obj,dt,varargin)
            % obj = CheckTimestep(obj,dt,varargin)
            % varargin(1) is desired CFL (< 1) or maximum iterations (> 1)
            % varargin(2) is desired dj_cutoff 
            %% Decimate mesh to achieve CFL condition for stability.
            % Takes a mesh and removes triangles and nodes to produce a mesh
            % that satisfies the given timestep requirements of the user by
            % incrementally modifying the triangulation nearby each edge
            % that violates the CFL.
            % kjr, chl, und, 2017
            % kjr, chl, und, 2018 <--updated for projected spaces
            % wjp, chl, und, 2019 <--updates to carry over slopes and
            % using the msh clean, and CFL calc type option
            
            type       = 0;       %<-- 0 == Haversine, 1 == CPP with correction factor
            desCFL     = 0.50;    %<-- set desired cfl (generally less than 0.80 for stable).
            desIt = inf;          %<-- desired number of iterations
            djc   = 0.25;
            if nargin > 2
                if varargin{1} > 1
                    desIt = varargin{1};
                else
                    desCFL = varargin{1};
                end
                if nargin == 4
                    djc = varargin{2};
                    %type = varargin{2}; 
                end
            end
            %%
            if ~isempty(obj.coord)
                % kjr 2018,10,17; Set up projected space imported from msh class
                global MAP_PROJECTION MAP_VAR_LIST MAP_COORDS
                MAP_PROJECTION = obj.proj ;
                MAP_VAR_LIST   = obj.mapvar ;
                MAP_COORDS     = obj.coord ;
            else
                lon_mi = min(obj.p(:,1)); lon_ma = max(obj.p(:,1));
                lat_mi = min(obj.p(:,2)); lat_ma = max(obj.p(:,2));
                m_proj('Trans','lon',[lon_mi lon_ma],'lat',[lat_mi lat_ma]) ;
            end
            
            %% Project..
            [obj.p(:,1),obj.p(:,2)] = m_ll2xy(obj.p(:,1),obj.p(:,2)); 
            
            %% Make bathy interpolant
            F = scatteredInterpolant(obj.p(:,1),obj.p(:,2),obj.b,...
                                     'linear','nearest');             
            if ~isempty(obj.bx)
                Fx = scatteredInterpolant(obj.p(:,1),obj.p(:,2),obj.bx,...
                                          'linear','nearest');
                Fy = scatteredInterpolant(obj.p(:,1),obj.p(:,2),obj.by,...
                                          'linear','nearest');
            end  
             
            % clear some things which cause error in renum
            obj.bx = []; obj.by = []; obj.f13 = [];
            obj.bd = []; obj.op = [];
            
            %% Delete CFL violations
            it  = 0;
            CFL = 999;
            if ~isempty(obj.pfix)
                [pf(:,1),pf(:,2)] = m_ll2xy(obj.pfix(:,1),obj.pfix(:,2));
            end
            con = 9;
			badnump = 1e10;
            tic
            while 1
                [obj.p(:,1),obj.p(:,2)] = m_xy2ll(obj.p(:,1),obj.p(:,2));
                CFL = CalcCFL(obj,dt,type);
                [obj.p(:,1),obj.p(:,2)] = m_ll2xy(obj.p(:,1),obj.p(:,2));
                bad = real(CFL) > desCFL;
                %if ~isempty(obj.pfix)
                %    fix = ismembertol(obj.p,pf,1e-6,'ByRows',true); 
                %    bad = bad & fix;
                %end
                display(['Number of CFL violations ',num2str(sum(bad))]);
                disp(['Max CFL is : ',num2str(max(real(CFL)))]);
                if it == desIt,   break; end
				badnum = sum(bad);
                if badnum == 0; break; end
                if badnump - badnum <= 0; con = con + 1; end
                it = it + 1;
                obj = DecimateTria(obj,bad);
                % Clean up the new mesh (already projected) without direct
                % smoothing (we have used local smooting in DecmimateTria
                obj.b = []; % (the bathy will cause error in renum)
                obj = clean(obj,'passive','proj',0,'pfix',pf,'con',con);
                obj.b = F(obj.p(:,1),obj.p(:,2));
            end
            toc
            disp(['Achieved max CFL of ',num2str(max(real(CFL))),...
                ' after ',num2str(it),' iterations.']);
            
            % add bathy back on
            obj.b = F(obj.p(:,1),obj.p(:,2));
            if exist('Fx','var')
                obj.bx = Fx(obj.p(:,1),obj.p(:,2));
                obj.by = Fy(obj.p(:,1),obj.p(:,2));
            end
            
            disp('Boundary and f13 info has been deleted...');
            disp('bathy and slope info is carried over');
            
            % find nans
            if ~isempty(find(isnan(obj.b), 1))
               warning('NaNs in bathy found')
            end
            
            % convert back to lat-lon wgs84
            [obj.p(:,1),obj.p(:,2)] = m_xy2ll(obj.p(:,1),obj.p(:,2));
            
            % Check Element order
            obj = CheckElementOrder(obj);
            return;
            
            function obj = DecimateTria(obj,bad)
                % form outer polygon of mesh for cleaning up.
                bnde = extdom_edges2(obj.t,obj.p);
                poly1 = extdom_polygon(bnde,obj.p,-1);
                poly_vec1 = cell2mat(poly1');
                [edges1]  = Get_poly_edges(poly_vec1);
                % delete the points that violate the cfl, locally retriangulate, and then locally smooth
                DTbase = delaunayTriangulation(obj.p(:,1),obj.p(:,2));
                DTbase.Points(find(bad),:) = [];
                pm   = DTbase.Points; tm = DTbase.ConnectivityList;
                pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/3;
                in1  = inpoly(pmid,poly_vec1,edges1);
                tm(~in1,:) = []; [pm,tm] = fixmesh(pm,tm);
                % Delete shitty boundary elements iteratively
                badbound = 1;
                while ~isempty(find(badbound, 1))
                    tq1 = gettrimeshquan(pm,tm);
                    % Get the elements that have a boundary bar
                    bnde = extdom_edges2(tm,pm);
                    bdnodes = unique(bnde(:));
                    vtoe = VertToEle(tm);
                    bele = unique(vtoe(:,bdnodes)); bele(bele == 0) = [];
                    tqbou = tq1.qm(bele); 
                    badbound = bele(tqbou < 0.1);
                    % Delete those boundary elements with quality < 0.1
                    tm(badbound,:) = [];
                end
                [pm,tm] = fixmesh(pm,tm);
                if sum(bad) < size(pm,1)*0.1
                    % find all the points nearby each "bad" point.
                    idx = ourKNNsearch(pm',obj.p(find(bad),:)',12);
                    idx = idx(:);
                    idx = unique(idx);
                    constr = setdiff((1:length(pm))',idx);
                    [pm,tm] = smoothmesh(pm,tm,constr,50,0.01);
                else
                    warning('Adapation would result in potentially catastrophic lose of connectivity, try adapting to a smaller timestep');
                end
                obj.p = pm; obj.t = tm;
            end
        end
        
        
        function [] = GatherStationsFromKML(obj,kmlfile,ofname,varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Read a kml file created in google maps to get coordinates of
            % stations that lie within the meshes extents.
            % (Note we can move the stations to reasonable locations off land etc in google
            % maps before exporting to kml).
            %
            %
            % Inputs:
            %     obj     - msh object
            %     kmlfile - kml filename dowloaded from William Pringle's tidal database.
            %     ofname  - filename of csv containing stations that fall within the mesh.
            %     const (optional)  - A cell-array of constitutents to validate against.
            %                         By default will be all major 8 (
            %                         {'M2','S2','N2','K2','K1','O1','P1','Q1'});
            %     Source (optional) - A cell-array of the source databases
            %                        in order of priority/reliability. By default it is:
            %                        sources = {'Truth_Pelagic','Truth_Shelf','NOAA','JMA','AUS_Tides',...
            %                        'KHAO','GESLA','UHSLC_Fast','Truth_Coast',...
            %                        'Various_Publications','ST727','IHO',};
            % Outputs:
            %
            % Author: William Pringle
            % Changes by Keith Roberts: merged code into msh class and added input parser for name/value pairs, 2018-03-22
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Variables to set
            if nargin < 3, error('Not enough input arguments!'); end
            
            if ~isempty(varargin)
                names = {'const','Source'};
                for ii = 1:length(names)
                    ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
                    if ~isempty(ind)
                        if ii == 1
                            const = varargin{ind*2};
                        else
                            Source = varargin{ind*2};
                        end
                    end
                end
                
            else
                % Set the constituents we wanna look at
                const = {'M2','S2','N2','K2','K1','O1','P1','Q1'};
                Source = {'Truth_Pelagic','Truth_Shelf','NOAA','JMA','AUS_Tides',...
                    'KHAO','GESLA','UHSLC_Fast','Truth_Coast',...
                    'Various_Publications','ST727','IHO',};
            end
            
            %% Reading the kml
            T = [];
            for kml = kmlfile
                kmlS = kml2struct(kml{1});
                % Getting rid of unnecessary shit
                T_n = struct2table(kmlS);
                T_n.Geometry = [];
                T_n.BoundingBox = [];
                %T.Description = [];
                T_n = T_n(:,[1 3 4 2 5]);
                T_n.Source = T_n.StyleUrl;
                T_n.StyleUrl = [];
                T = [T; T_n];
            end
            
            %% Convert description to amp and phases
            obs = NaN(height(T),length(const)*2);
            for i = 1:height(T)
                if strcmp(T.Source{i}(end-5:end),'0288D1') % light blue
                    T.Source{i} = 'Truth_Pelagic';
                elseif strcmp(T.Source{i}(end-5:end),'FFEA00') % yellow
                    T.Source{i} = 'Truth_Shelf';
                elseif strcmp(T.Source{i}(end-5:end),'0F9D58') % green
                    T.Source{i} = 'Truth_Coast';
                elseif strcmp(T.Source{i}(end-5:end),'757575') % grey
                    T.Source{i} = 'NOAA';
                elseif strcmp(T.Source{i}(end-5:end),'673AB7') % dark purple
                    T.Source{i} = 'JMA';
                elseif strcmp(T.Source{i}(end-5:end),'FF5252') % coral
                    T.Source{i} = 'KHAO';
                elseif strcmp(T.Source{i}(end-5:end),'7CB342') % light green
                    T.Source{i} = 'GESLA';
                elseif strcmp(T.Source{i}(end-5:end),'A52714') % medium red
                    T.Source{i} = 'UHSLC_Fast';
                elseif strcmp(T.Source{i}(end-5:end),'FFD600') % austrlian gold
                    T.Source{i} = 'AUS_Tides';
                elseif strcmp(T.Source{i}(end-5:end),'880E4F') % dark blue
                    T.Source{i} = 'ST727';
                elseif strcmp(T.Source{i}(end-5:end),'3949AB') % dark blue
                    T.Source{i} = 'IHO';
                elseif strcmp(T.Source{i}(end-5:end),'9C27B0') % light purple
                    T.Source{i} = 'Various_Publications';
                end
                Descrip = T.Description(i); Descrip = Descrip{1};
                cc = 0;
                for c = const
                    cc = cc + 1;
                    % Get position in description of the constituent
                    Positions = regexp(Descrip,[c{1} '_']);
                    Positions = find(~cellfun(@isempty,Positions));
                    if isempty(Positions); continue; end
                    if length(Positions) > 2
                        % M2 seasonal shit
                        Positions(1:2) = [];
                    end
                    amp = Descrip(Positions(1));
                    s = strfind(amp,' ');
                    obs(i,cc*2-1) = str2double(amp{1}(s{1}+1:end));
                    phs = Descrip(Positions(2));
                    s = strfind(phs,' ');
                    obs(i,cc*2) = str2double(phs{1}(s{1}+1:end));
                end
            end
            % Make names
            for ii = 1:length(const)
                const_names{2*ii-1} = [const{ii} '_amp'];
                const_names{2*ii} = [const{ii} '_phs'];
            end
            T1 = array2table(obs,'VariableNames',const_names);
            T.Description = [];
            % Combine
            T = [T T1];
            
            %% Only keep stations within our domain
            bnde = extdom_edges2(obj.t,obj.p);
            poly = extdom_polygon(bnde,obj.p,0);
            
            k=0; poly_vec=[];
            for i = 1 : length(poly)
                for ii = 1 : length(poly{i})
                    k = k + 1;
                    poly_vec(k,:) = poly{i}(ii,:);
                end
                k = k + 1;
                poly_vec(k,:) = [NaN,NaN];
            end
            edges = Get_poly_edges(poly_vec);
            in = inpoly([T.Lon,T.Lat],poly_vec,edges);
            T = T(in,:);
            %% Delete duplicates
            radius = 0.08/60; % only take one station within 1 min box
            [~,IA] = uniquetol([T.Lon T.Lat],radius,'ByRows',true,...
                'DataScale',1,'OutputAllIndices',true);
            % Just initialise T_new randomly
            T_new = T(1:length(IA),:);
            delete = []; nn = 0;
            for ii = 1:length(IA)
                T_temp = T(IA{ii},:);
                % We want to pick the source in heirachial order
                for s = Source
                    exp_s = regexp(T_temp.Source,s{1});
                    exp_s = find(~cellfun(@isempty,exp_s));
                    if ~isempty(exp_s)
                        if length(exp_s) > 1
                        end
                        T_new(ii,:) = T_temp(exp_s,:);
                        break
                    end
                end
                if isempty(exp_s)
                    nn = nn + 1;
                    delete(nn) = ii;
                end
            end
            T_new(delete,:) = [];
            % sorting
            [~,IA] = sort(T_new.Name);
            T_new = T_new(IA,:);
            [~,IA] = sort(T_new.Source);
            T_new = T_new(IA,:);
            
            %% Write to the .csv
            % Output table of the names and positions
            writetable(T_new(IA,:),ofname)
            %             %% Populate the f15
            %             % elevation
            %             obj.f15.nstae = obj.f15.nstae + numel(find(Sta_type(:,1)));
            %             obj.f15.elvstaloc = [obj.f15.elvstaloc;
            %                 [Sta_lon(Sta_type(:,1) == 1) ...
            %                 Sta_lat(Sta_type(:,1) == 1)]];
            %             obj.f15.elvstaname = [obj.f15.elvstaname;
            %                 strcat(Sta_name(Sta_type(:,1) == 1),' ID:',...
            %                 Sta_ID(Sta_type(:,1) == 1))];
            
        end
        
        function [ef,efx,efy]=reconstructEdgefx(obj,efdx)
            % Given a msh object, reconstruct the edge function resolution
            % with a resolution equal to efdx in WGS84 degrees.
            efdx = efdx/111e3;
            
            TR = triangulation(obj.t,obj.p(:,1),obj.p(:,2));
            [~,cr] = circumcenter(TR);
            
            for i = 1  : length(obj.t)
                cl = cr(i) ;
                for j = 1 :  3
                    z(obj.t(i,j)) = cl;
                end
            end
            
            bbox = [min(obj.p); max(obj.p)]';
            [efx,efy] = ndgrid(bbox(1,1):efdx:bbox(1,2),...
                bbox(2,1):efdx:bbox(2,2));
            F = scatteredInterpolant(obj.p(:,1),obj.p(:,2),z','linear','linear');
            
            ef = F(efx,efy);
            
            bnde = extdom_edges2(obj.t,obj.p);
            poly1=extdom_polygon(bnde,obj.p,1);
            for i = 1 : length(poly1)
                poly1{i} = [poly1{i} ; NaN NaN];
            end
            poly_vec1=cell2mat(poly1');
            edges1=Get_poly_edges(poly_vec1);
            in = inpoly([efx(:),efy(:)],poly_vec1,edges1);
            ef(~in) = NaN;
        end
        
        function [centroids]=baryc(obj)
            centroids = (obj.p(obj.t(:,1),:)+obj.p(obj.t(:,2),:)+obj.p(obj.t(:,3),:))/3;
        end
        
        function [obj3] = minus(obj1,obj2)
            % Given two overlapping triangulations described in msh objects
            % obj1 and obj2, remove the triangles in obj1 with centroids
            % that are in the boundary of obj2.
            % Here the notion of "in" is defined as being with the polygon
            % that is the boundary of obj2.
            %
            % NOTE: The first mesh must define a region that inclues the region
            % of obj2 but not the entitery, otherwise the algorithm removes
            % all of obj1 and segfaults.
            % Author: keith, und, chl, 2018
            bnde = extdom_edges2(obj2.t,obj2.p);
            poly = extdom_polygon(bnde,obj2.p,-1);
            polyc= cell2mat(poly');
            edges=Get_poly_edges(polyc);
            
            centroids = baryc(obj1) ;
            in = inpoly(centroids,polyc,edges);
            
            tempp = obj1.p;
            tempt = obj1.t;
            
            tempt(in,:) = [] ;
            [tempp,tempt]=fixmesh(tempp,tempt) ;
            
            obj3 = msh();
            obj3.p = tempp; obj3.t = tempt;
            obj3 = renum(obj3);
            
        end
        
        
        function obj = bound_con_int(obj,nvl) %nvu,
            %  bound_con_int updates all interior nodes that are connected
            %  to more than 7 nodes and less than 21 so that every node is
            %  connected to at most 7 nodes. It also makes sure there are no
            %  small connectivity less than 5 nodes. A springing routine is used locally
            %  to help smooth out the updates.  Note that specifying NVL one can
            %  fix nodes with connectivity of NVL or higher and leave the others.
            %
            %  NOTE:  fem_struct.nei should already be a component, if not one is
            %         computed.
            %
            %
            %  Variables
            %  fem_struct -- is the finite element structure from the opnml suite.
            %  fem -- the updated finite element structure.
            %        Note:  fem.x,fem.y,fem.z,fem.e, fem.nei, and fem.ar
            %                  get updated.
            %               fem.bnd does not need to be updated.
            %  nvl -- lower limit of connectivity to fix should be 8 <= nvl <= 20
            %  nvu -- upper limit of connectivity to fix should be 4 usually
            %
            %  Usage -- fem = reduce_con_int(fem_struct);
            %
            %  Name: reduce_con_int.m
            %  Written by: Ben Holladay (SEAP Student 2004)
            %  Date: June 22,2004
            %  Modified:  Aug. 30, 2004, Chris Massey
            %             Jan. 10, 2006, Chris Massey
            %             Aug. 15, 2006, Chris Massey  - Cleaned up comments.
            %             Sept. 19, 2006, Chris Massey -- fixed bug
            %             Aug. 22, 2007, Chris Massey -- fixed bug
            %             Apr. 25, 2018, Chris Massey -- Added nvl option to fix nodes
            %                with connectivity higher than or equal to nvl.
            %             Apr. 26, 2018, Keith Roberts--Merged into msh class for
            %             OceanMesh software
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % convert msh obj to fem_struct
            bnde=extdom_edges2(obj.t,obj.p);
            fem_struct.x = obj.p(:,1);
            fem_struct.y = obj.p(:,2);
            fem_struct.e = obj.t;
            if ~isempty(obj.b)
                fem_struct.z = obj.b;
            else
                fem_struct.z = obj.p(:,1)*0;
            end
            fem_struct.bnd = bnde;
            fem_struct.ar  = fem_struct.z*0.0;
            fem_struct.name= 'temp';
            if ~isempty(obj.f15) || ~isempty(obj.f13)
                warning('Reduce nodal connectivity will erase your nodal attributes and control file. Go on?');
                pause
            end
            
            %Ensure the appropriate number of input arguements is given, checks that the
            %fem_struct is valid.
            if nargin == 0
                error('Not enough input arguments; need a valid fem_struct.');
            end
            if ~is_valid_struct(fem_struct)
                error('Input argument to reduce_con_bnd must be a valid fem_struct.');
            end
            %             % Make sure that the upper limit of connectivity count falls
            %             % no larger than 4 nodes
            %             if nargin < 2
            %                 if nvl > 4
            %                     disp('Connectivity Upper Limit must not be greater than 4. Resetting');
            %                     nvu = 4;
            %                 end
            %             else
            %                 nvu = 4;
            %             end
            
            %Make sure that the lower limit of connecitivity count falls between
            % 8 and 20 nodes
            if nargin == 2
                if nvl < 8
                    disp('Connectivity Lower Limit must be at least 8.  Resetting');
                    nvl = 8;
                end
                if nvl > 20
                    disp('Connectivity Lower Limit must not be greater than 20. Resetting');
                    nvl = 20;
                end
            else
                nvl = 8;
            end
            
            %Sets the intial fem_struct variables and determines the max
            %connectivity. Displays errors if the connectivity is to high.
            
            % List of all nodes not on the boundary;
            tmpnodes = setdiff(1:1:length(fem_struct.x),unique(fem_struct.bnd(:)));
            bndrynodes = unique(fem_struct.bnd(:));
            
            nelems = size(fem_struct.e,1);
            
            % Determine if fem_struct.nei is present, if not added it.
            try
                nc = size(fem_struct.nei(tmpnodes,:),2);
            catch
                disp(' ');
                disp('A neighbor list was not present in the finite element structure.');
                disp('One is being added now.');
                fem_struct.nei = ele2nei(fem_struct.e,fem_struct.x,fem_struct.y);
                disp(' ');
                disp('The neigbor list was successfully added.');
                nc = size(fem_struct.nei(tmpnodes,:),2);
            end
            
            nelems_orig = nelems;
            nnodes_orig = length(fem_struct.x);
            nc_orig = nc;
            
            [~,J1] = find(fem_struct.nei(bndrynodes,:)~=0);
            highbndry = max(J1);
            clear J1
            
            %nsmall = any(fem_struct.nei(j,nvu+1));
            
            if nc > 21
                error(['The connectivity is too high and must be brought down to at least',...
                    '20 by hand.']);
            end
            
            if nc <= nvl-1
                disp(['The grids nodal connectivity is not higher than ',num2str(nvl),', nothing will be done.']);
                disp('Returning the original mesh.');
                return
            end
            
            disp(' ')
            disp('NOTE: NO UPDATING IS DONE FOR BOUNDARY NODES.');
            disp(' ')
            
            %Begins loop to update nodes.
            tempr = size(fem_struct.nei,2);
            if tempr < nvl-1 %7
                tempr = 0;
            else
                tmpr = length(find(fem_struct.nei(tmpnodes,nvl)~=0));
            end
            imax = 20;i = 1;
            i2flag = 1;
            fixedint=[];
            disp('Beginning the update loop')
            disp(' ')
            while tempr ~= 0 && i <= imax
                disp(['Pass number ',num2str(i,'%4.0f'),' through the loop'])
                
                %tmpnodes2 = tmpnodes;
                %iupdnn=find(fem_struct.nei(tmpnodes,nvl)~=0);  %nodes with high connectivity
                %tmpnodes = tmpnodes(iupdnn);
                %clear tmpnodes2
                nnodes = size(fem_struct.nei(tmpnodes,:),1);
                
                j2 = 1;
                while j2 <= nnodes
                    jj = tmpnodes(j2);
                    temp2 = fem_struct.nei(jj,nvl); %8;
                    if  temp2 ~= 0
                        testnei = fem_struct.nei(jj,:);
                        % TCM -- Begin -- 09/19/2006
                        %temp1 = sum(ismember((fem_struct.nei(j,:)),0)); TCM
                        %tempnc = nc - temp1; TCM
                        tempnc = sum(~ismember((fem_struct.nei(jj,:)),0)); %TCM
                        % TCM -- END -- 09/19/2006
                        
                        disp(['Updating node number ',num2str(jj),' using',' fix ',num2str(tempnc)]);
                        switch tempnc
                            case {13,14,15,16,17,18,19,20}
                                fem_struct = update1320nbr(fem_struct,jj);
                            case 12
                                fem_struct = update12nbr(fem_struct,jj);
                            case 11
                                fem_struct = update11nbr(fem_struct,jj);
                            case 10
                                fem_struct = update10nbr(fem_struct,jj);
                            case 9
                                fem_struct = update9nbr(fem_struct,jj);
                            case 8
                                fem_struct = update8nbr(fem_struct,jj);
                            otherwise
                                disp('Already fixed.')
                        end
                        if sum(testnei) == sum(fem_struct.nei(jj,:)) %?? What is this test?
                            fixedint(i2flag) = jj;
                            temp1 = find(fem_struct.nei(jj,:) ~= 0,1,'last');
                            fixedintnei(i2flag) = temp1;
                            i2flag = i2flag + 1;
                        end
                    end
                    if size(fem_struct.nei,2) < nvl %8
                        j2 = nnodes + 1;
                    else
                        j2 = j2 + 1;
                    end
                end
                
                %Check if any nodes where added to the do not update lists.
                try
                    fixedint = fixedint;
                    fixedintnei = fixedintnei;
                catch
                    fixedint = 1:0;
                    fixedintnei = 1:0;
                end
                
                %Create do not update list for interior nodes.
                
                %newnei = fem_struct.nei(fixedint,:);
                %temp1 = ismember(fem_struct.nei(fixedint,:),0);
                %temp2 = sum(temp1,2);
                %temp3 = [(size(fem_struct.nei,2))-temp2];
                %temp4 = fixedintnei' - temp3;
                %temp5 = find(temp4 < 0);
                %fixedint(temp5) = [];
                %fixedintnei(temp5) = [];
                %temp6 = size(temp5,1);
                %i2flag = i2flag - temp6;
                
                int = setdiff(1:length(fem_struct.x),bndrynodes);
                
                tmpnodes = setdiff(int,fixedint);
                
                i = i + 1;
                %tmpnodes = setdiff(1:1:length(fem_struct.x),unique(fem_struct.bnd(:)));
                if size(fem_struct.nei,2) > nvl-1
                    tempr = length(find(fem_struct.nei(tmpnodes,nvl) ~= 0));
                else
                    tempr = 0;
                end
            end
            
            
            
            %Creates the output structure.
            fem = fem_struct;
            [~,J1] = find(fem_struct.nei(tmpnodes,:)~=0);
            highint = max(J1);
            
            fem.nei = fem.nei(:,1:max(highint,highbndry));
            fem = el_areas(fem);
            
            
            % Display a message that we reached the maximum number of loop iterations
            % prior to minimizing the interior nodal connectivity list
            if i>imax && highint > nvl-1 %7
                disp(' ')
                disp(['Maximum loop iteration count of ',num2str(imax),' was',...
                    ' exceeded before all interior nodes']);
                disp('reached their miniminal connecitvity.');
            end
            
            % Display a summary of mesh updates.
            disp(' ')
            disp([num2str(size(fem.e,1)-nelems_orig),' elements have been added']);
            disp([num2str(length(fem.x)-nnodes_orig),' nodes have been added']);
            disp(['Maximum interior nodal connectivity was reduced from ',...
                num2str(nc_orig), ' to ',num2str(highint),'.']);
            disp(['The maximum nodal connectivity for a boundary node',...
                ' is ',num2str(highbndry),'.']);
            disp(' ');
            
            
            % kjr convert back to msh obj.
            obj.p = [fem_struct.x,fem_struct.y];
            if ~isempty(obj.b)
                obj.b = fem_struct.z;
            end
            obj.t = fem_struct.e;
        end
        
        function obj = flipEdge(obj,j,k)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Line_swap is mesh refining tool that swaps an edge that two elements
            % share. It reconnects the elements so that the edge they share has been
            % flipped. This routine only operates on one line at a time and can only be
            % used if the two elements share an edge. This routine can operate on any
            % elemnet, even boundry elements, unless the new elements formed by the line
            % swap would overlap existing elements. No nodes are moved by this routine.
            % The fem.e and fem.ar are the only fields that are updated. If no .nei file
            % is present then none will be generated; however, if one is present then
            % it will be updated to reflect the new connectivity.
            %
            % Calls: is_valid_struct.m
            %
            % Usage: fem = line_swap(fem_struct,j,k);
            %
            % Variables:
            %  fem -- the new split finite element mesh.
            %  fem_struct -- the finite element grid structure from the opnml suite.
            %          Note: fem_struct.nei will not be generated if not present
            %  j -- element number for one element used in the line swap.
            %  k -- element number for one element used in the line swap.
            %          Note: j and k can be interchanged.
            %
            % Filename: line_swap.m
            % Created by: Ben Holladay
            % Date: June 8, 2005
            % Last Modified:
            %       Oct. 26, 2007 -- Chris Massey, NRL 7322
            %       Changed when line swap can not be performed due to mesh tangling,
            %       from an error call to a display message and return original fem.
            % Last Modified:
            %       Added into msh class for OceanMesh2D software--kjr, chl, und,
            %       April 2018.
            %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % convert msh obj to fem_struct
            bnde=extdom_edges2(obj.t,obj.p);
            fem_struct.x = obj.p(:,1);
            fem_struct.y = obj.p(:,2);
            fem_struct.e = obj.t;
            if ~isempty(obj.b)
                fem_struct.z = obj.b;
            else
                fem_struct.z = obj.p(:,1)*0;
            end
            fem_struct.bnd = bnde;
            fem_struct.ar  = obj.b*0.0;
            fem_struct.name= 'temp';
            if ~isempty(obj.f15) || ~isempty(obj.f13)
                warning('Reduce nodal connectivity will erase your nodal attributes and control file. Go on?');
                pause
            end
            
            %Ensure the apropriate number of input arguements is given, checks that the
            %fem_struct is valid.
            if nargin == 0
                error('Not enough input arguments; need a fem_struct and two element numbers.');
            end
            if ~is_valid_struct(fem_struct)
                error('Input argument to line swap must be a valid fem_struct.');
            end
            if nargin == 1
                error(['Input arguments must include a valid fem_struct and two element ' , ...
                    'numbers that share an edge.']);
            elseif nargin == 2
                error('Two element numbers must be given.');
            elseif nargin == 3
                if j <= 0 | j > size(fem_struct.e,1) | j ~= floor(j)
                    error(['Element numbers must positive integers that are less than ',...
                        'the maximum number of elements.']);
                end
                if k <= 0 | k > size(fem_struct.e,1) | k ~= floor(k)
                    error(['Element numbers must positive integers that are less than ',...
                        'the maximum number of elements.']);
                end
            else
                error('Too many input arguments.');
            end
            
            
            %Sets the intial fem_struct variables.
            x = fem_struct.x;
            y = fem_struct.y;
            enodes = fem_struct.e;
            ar = fem_struct.ar;
            
            %Determine if the elements share a common edge.
            elems = enodes([j;k],:);
            comp = unique(intersect(elems(1,:),elems(2,:)));
            if length(comp) < 2
                error('The elements do not share a common edge.');
            end
            
            %Indentify each node for swapping.
            [nodes,ti,tj] = unique(elems);
            temp1 = setdiff((1:6),ti);
            temp1 = elems(temp1);
            temp2 = setdiff(nodes,temp1);
            
            %Determine if the angle between the elements is too big for line_swap to
            %effectively operate on it.
            a2 = (x(enodes([j;k],3))-x(enodes([j;k],2))).^2+(y(enodes([j;k],3))-y(enodes([j;k],2))).^2;
            b2 = (x(enodes([j;k],1))-x(enodes([j;k],3))).^2+(y(enodes([j;k],1))-y(enodes([j;k],3))).^2;
            c2 = (x(enodes([j;k],2))-x(enodes([j;k],1))).^2+(y(enodes([j;k],2))-y(enodes([j;k],1))).^2;
            A = (180/pi)*acos((b2+c2-a2)./(2*sqrt(b2).*sqrt(c2)));
            B = (180/pi)*acos((c2+a2-b2)./(2*sqrt(c2).*sqrt(a2)));
            C = (180/pi)*acos((a2+b2-c2)./(2*sqrt(a2).*sqrt(b2)));
            ang = [A,B,C];
            temp3 = find(temp1(1) == elems);
            temp4 = find(temp1(2) == elems);
            temp5 = ang(temp3);
            temp6 = ang(temp4);
            temp5 = sum(temp5);
            temp6 = sum(temp6);
            if temp5 >= 180 | temp6 >= 180
                disp(['The line cannot be swapped because the new elements would '...
                    'create overlapping elements.']);
                fem = fem_struct;
                return
            elseif temp5 >= 135 | temp6 >= 135
                disp(['Warning: The elements created will contain large angles.']);
            end
            clear a2 b2 c2 A B C ang
            
            %Swap the lines in the enodes list.
            temp3 = setdiff(elems(1,:),temp1);
            temp4 = setdiff(elems(2,:),temp1);
            temp5 = find(elems(1,:) == temp1(1));
            elems(1,temp5) = temp4;
            temp6 = find(elems(2,:) == temp1(2));
            elems(2,temp6) = temp3;
            clear temp5 temp6;
            
            %Recomputes the areas for the new elements.
            xnodes = x(elems);
            ynodes = y(elems);
            temparea = 0.5*(xnodes(:,1).*(ynodes(:,2)-ynodes(:,3))+xnodes(:,2).*...
                (ynodes(:,3)-ynodes(:,1))+xnodes(:,3).*(ynodes(:,1)-ynodes(:,2)));
            clear xnodes ynodes x y;
            
            %Addes the new elements and areas to the global lists.
            enodes([j;k,],:) = elems;
            ar([j;k]) = temparea(:);
            clear elems temparea;
            
            %Create the output structure.
            fem = fem_struct;
            fem.e = enodes;
            fem.ar = ar;
            
            %Updates the nei if present.
            try
                nei = fem_struct.nei;
                nei = [nei,zeros((size(nei,1)),1)];
                tempnei = nei([temp1(:);temp2(:)],:);
                temp3 = find(tempnei(1,:) == temp1(2));
                temp4 = setdiff(1:(size(tempnei,2)),temp3);
                tempnei(1,:) = ([tempnei(1,temp4),0]);
                temp3 = find(tempnei(2,:) == temp1(1));
                temp4 = setdiff(1:(size(tempnei,2)),temp3);
                tempnei(2,:) = ([tempnei(2,temp4),0]);
                temp3 = find(tempnei(3,:) == temp1(1) | tempnei(3,:) == temp1(2));
                if abs(temp3(1)-temp3(2)) == 1
                    tempnei(3,:) = [tempnei(3,(1:temp3(1))),temp2(2),tempnei(3,((temp3(2)):(end-1)))];
                else
                    tempnei(3,:) = [temp2(2),tempnei(3,(1:end-1))];
                end
                temp3 = find(tempnei(4,:) == temp1(1) | tempnei(4,:) == temp1(2));
                if abs(temp3(1)-temp3(2)) == 1
                    tempnei(4,:) = [tempnei(4,(1:temp3(1))),temp2(1),tempnei(4,((temp3(2)):(end-1)))];
                else
                    tempnei(4,:) = [temp2(1),tempnei(4,(1:end-1))];
                end
                nei([temp1(:);temp2(:)],:) = tempnei;
                temp3 = find(nei(:,end) ~= 0);
                if isempty(temp3) == 1
                    nei = nei(:,(1:end-1));
                end
                temp3 = find(nei(:,end) ~= 0);
                if isempty(temp3) == 1
                    nei = nei(:,(1:end-1));
                end
                fem.nei = nei;
            catch
            end
            
            % kjr convert back to msh obj.
            obj.p = []; obj.b = []; obj.t =[];
            obj.p = [fem_struct.x,fem_struct.y];
            obj.b = fem_struct.z;
            obj.t = fem_struct.e;
        end
        
        function mfixed = MergeFP(muw,mfp) 
            %%%%%%%
            % Merges a watertight mesh muw with a mesh that contains both 
            % over and underwater sections mfp. This method assumes mfp has
            % been created by edge locking the shoreline and using fixed
            % points. 
            %
            % INPUTS: 
            % muw: watertight msh object 
            % mfp: watertight msh with the floodplain using edge locking
            %
            % OUTPUTS: 
            % mfixed: a final merged msh that seamlessly mates with the
            % underwater mesh.
            % kjr, April 2019
            
            m3 = mfp - muw ;
                        
            dt = delaunayTriangulation(muw.p) ;
            
            dt.Points(end+(1:length(m3.p)),:) = m3.p ;
            
            tmp = mfp;
            
            tmp.p = dt.Points ; tmp.t = dt.ConnectivityList; 
            
            % delete points outside mfp and muw 
            bnde = extdom_edges2(mfp.t,mfp.p) ;
            polyt=extdom_polygon(bnde,mfp.p,-1,0,20) ;
            polyt=cell2mat(polyt') ;
            
            ee   = Get_poly_edges(polyt) ;
            bc   = baryc(tmp) ;
            infp   = inpoly(bc,polyt,ee) ;
            
            bnde = extdom_edges2(muw.t,muw.p) ;
            polyt=extdom_polygon(bnde,muw.p,-1,0,20) ;
            polyt=cell2mat(polyt') ;
            
            ee   = Get_poly_edges(polyt) ;
            bc   = baryc(tmp) ;
            inuw   = inpoly(bc,polyt,ee) ;
                        
            tmp.t(~infp & ~inuw,:) = [] ; 
                       
            [tmp.p,tmp.t] = fixmesh(tmp.p,tmp.t) ; 
            
            mfixed = tmp ;
            
            mfixed = renum(mfixed) ; 
            
            
        end
        
        function [pfix,egfix] = extractFixedConstraints(obj)
            %%%%%%%
            % Extract boundary of mesh in no order.
            % INPUTS: msh_obj 
            % OUTPUTS: the points and edges of the mesh 
            % NOTE: You can visualize these constraints by 
            % drawedge2(pfix,egfix); 
            % kjr, April 2019
            [egfix,pfix] = extdom_edges2(obj.t,obj.p) ;
            if ~isempty(obj.op) && obj.op.nope > 0
                disp('detected ocean boundary, removing these fixed points')
                ocean = unique(obj.op.nbdv(:));
                [~,IA] = intersect(egfix(:,1),ocean);
                [~,IB] = intersect(egfix(:,2),ocean);
                IC = unique([IA;IB]);
                egfix(IC,:) = []; 
                pfix = obj.p(unique(egfix(:)),:); 
            end
            egfix = renumberEdges(egfix) ;
        end
        
        function boundary = getBoundaryOfMesh(obj) 
            %%%%%%%
            % Returns the boundary of the mesh in a winding order as a
            % NaN-delimited vector 
            % INPUTS: msh_obj 
            % OUTPUTS: NaN-delimited vector in a "walking" order.
            % kjr, April 2019
            bnde = extdom_edges2(obj.t,obj.p) ; 
            try
                boundary = extdom_polygon(bnde,obj.p,-1) ;
            catch
                warning('ALERT: Boundary of mesh is not walkable. Returning polylines.');
                boundary = extdom_polygon(bnde,obj.p,-1,1) ;
            end
            boundary = cell2mat(boundary'); 
        end
        
        function obj = pruneOverlandMesh(obj,elev,djc) 
            %%%%%%%
            % Removes overland extent greater than certain elevation about 
            % geoid of the mesh.
            % INPUTS: msh_obj and height above the geoid whereby elements with
            % an average depth greater than elev are pruned from the mesh. 
            % OUTPUS: a msh_obj with the overland extents removed. 
            % april 2019, kjr
            if nargin < 3
               djc = 0.25; 
            end
            [obj.p(:,1),obj.p(:,2)] =  m_ll2xy(obj.p(:,1),obj.p(:,2)); 
            Fb = scatteredInterpolant(obj.p(:,1),obj.p(:,2),obj.b,'linear','nearest') ; 
            c = (obj.b(obj.t(:,1),:)+obj.b(obj.t(:,2),:)+obj.b(obj.t(:,3),:))/3;
            obj.t(c < -elev,:) = []; 
            [pp,tt] = fixmesh(obj.p,obj.t) ; 
            obj.p = pp; obj.t = tt; 
            obj = Make_Mesh_Boundaries_Traversable(obj,djc,1);
            obj = renum(obj) ;
            obj.b = Fb(obj.p) ; 
            [obj.p(:,1),obj.p(:,2)] =  m_xy2ll(obj.p(:,1),obj.p(:,2)); 
        end
        
        function obj = interpFP(obj,gdat,muw,gdatuw,minb,CAN)
            %%%%%%%
            % Interpolate topography onto a mesh with floodplain
            % Interpolates bathymetry underwater using a grid-scale
            % averaging with a minimum depth of 1-m. 
            % Then it uses a 3x multipler for grid-scale averaging
            % to interpolate overland topography. 
            % Finally, it updates the bathymetry on the outputted msh
            % to match the underwater and the overland sections. 
            %
            %%%%%%%
            % INPUTS
            % obj: msh_obj to interpolate topograhy/bathy on
            % gdat:geodata obj that contains DEM (could be cell-array of
            % gdats)
            % muw: msh_obj of only the watertight portion of the domain
            %
            %%%%%%%
            % OUTPUTS: 
            % a msh_obj with bathy/topo on its vertices 
            %
            % kjr, april 2019

            % parsing some inputs (or set to default)
            if nargin < 5
                minb = 1;
            end
            if nargin < 6
                CAN = 3;
            end
               
            bnde = extdom_edges2(muw.t,muw.p) ; 
            bou = extdom_polygon(bnde,muw.p,-1) ; 
            bou = cell2mat(bou') ; 
            
            dmy1 = obj; % land + uw
            
            ee = Get_poly_edges(bou); 
            in = inpoly(dmy1.p,bou,ee) ; 
            
            % get the "on"
            on = ismembertol(dmy1.p,muw.p,1e-5,'ByRows',true);
            
            obj.b = (1:length(obj.p(:,1)))'*0 ; 
            
            if ~iscell(gdat); gdat = {gdat}; end
            
            for i = 1 : length(gdat)
                if isempty(gdat{i}.Fb) 
                   disp(['Entry ',num2str(i), ' does not contain DEM, skipping']);
                   continue 
                end
                if i == 1
                    in2 = true(size(obj.b));
                else
                    in2 = inpoly(dmy1.p,gdat{i}.boubox(1:end-1,:)) ;
                end
                dmy1 = interp(obj,gdatuw(i),'type','depth','ignoreOL',1); 
                
                dmy2 = interp(obj,gdat(i),'type','depth','N',CAN) ; % use smooth overland
                
                if ~isempty(gdat{i}.mainlandb) || ~isempty(gdat{i}.innerb)
                    riverbound = [];
                    if ~isempty(gdat{i}.mainlandb)
                        notlakem = find(~contains(gdat{i}.mainlandb_type,'lake'));
                        isnan1 = find(isnan(gdat{i}.mainlandb(:,1)));    
                        for l = 1:length(notlakem)
                            if notlakem(l) == 1; isn1 = 1; else
                            isn1 = isnan1(notlakem(l)-1)+1; end
                            isn2 = isnan1(notlakem(l))-1;
                            riverbound = [riverbound; ...
                                          gdat{i}.mainlandb(isn1:isn2,:)];
                        end
                    end
                    if ~isempty(isempty(gdat{i}.innerb))
                        notlakei = find(~contains(gdat{i}.innerb_type,'lake'));
                        isnan1 = find(isnan(gdat{i}.innerb(:,1)));    
                        for l = 1:length(notlakei)
                            if notlakei(l) == 1; isn1 = 1; else
                            isn1 = isnan1(notlakei(l)-1)+1; end
                            isn2 = isnan1(notlakei(l))-1;
                            riverbound = [riverbound; ...
                                          gdat{i}.innerb(isn1:isn2,:)];
                        end
                    end

                    % kjr 2018,10,17; Set up projected space imported from msh class
                    dmyriver = dmy2; 
                    dmyriver.p = [dmyriver.p; riverbound(:,1:2)];
                    setProj(dmyriver,1,obj.proj.name)         
                    [rvb(:,1),rvb(:,2)] = ...
                                  m_ll2xy(riverbound(:,1),riverbound(:,2));
                    [dmp(:,1),dmp(:,2)] = m_ll2xy(dmy2.p(:,1),dmy2.p(:,2));                                    
                    F = scatteredInterpolant(rvb(:,1),rvb(:,2),...
                                      riverbound(:,3),'natural','nearest');
                    offset = F(dmp);
                    % change above land
                    dmy2.b = min(0,dmy2.b + offset); 
                    % change minb based on river height
                    minb = max(1,dmy1.b*0 + minb - offset);
                end
                uw = in | on | dmy1.b > 0;
                dmy1.b = max(dmy1.b,minb); % bound the depth by minb
                obj.b(uw & in2,1) = dmy1.b(uw & in2,1) ;
                obj.b(~uw & in2,1) = dmy2.b(~uw & in2,1) ;
                
            end
          
        end
        
        
    end % end methods
    
    
    
end % end class
