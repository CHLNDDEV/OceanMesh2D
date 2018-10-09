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
    end
    
    methods
        % constructor/read mesh into class.
        function obj = msh(fname,type)
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
        function h = plot(obj,type,proj,projection,bou)
            if nargin < 3
                proj = 1;
            end
            if nargin == 5
               if numel(bou) == 4
                    % i.e. is a bounding box
                    bou = [bou(1,1) bou(2,1);
                           bou(1,1) bou(2,2); ...
                           bou(1,2) bou(2,2);
                           bou(1,2) bou(2,1); ...
                           bou(1,1) bou(2,1)]; 
               end
               % Get a subset given by bou
               obj = ExtractSubDomain(obj,bou);
            end
            if proj
                if nargin < 4 || isempty(projection)
                    projection = 'Transverse Mercator';
                end
                if contains(projection,'ster')
                    m_proj(projection,'lon',mean(obj.p(:,1)),...
                           'lat',mean(obj.p(:,2)),...
                           'radius',0.6*max(max(obj.p)-min(obj.p)))  ;                    
                else
                    m_proj(projection,...
                           'lon',[min(obj.p(:,1)),max(obj.p(:,1))],...
                           'lat',[min(obj.p(:,2)),max(obj.p(:,2))])  ;
                end
            end
            logaxis = 0; numticks = 10;
            if strcmp(type(max(1,end-2):end),'log')
                logaxis = 1; type = type(1:end-3);
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
                                    obj.p(obj.op.nbdv(1:obj.op.nvdll(nb),nb),2),'b-','linewi',1.2);
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
                        q = obj.b;
                    end 
%                     if nargin == 5
%                         % Trick to full it out white space
%                         lon = linspace(min(bou(:,1)),max(bou(:,1)),500);
%                         lat = linspace(min(bou(:,2)),max(bou(:,2)),500);
%                         [lon,lat] = meshgrid(lon,lat);
%                         F = scatteredInterpolant(obj.p(:,1),obj.p(:,2),...
%                                                  q,'linear','nearest');
%                         V = F(lon,lat); V(25:end-25,25:end-25) = NaN;
%                         if proj
%                            m_pcolor(lon,lat,V); shading interp
%                         else
%                            pcolor(lon,lat,V); shading interp
%                         end
%                     end 
                    if proj
                        m_trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                        %m_trisurf(obj.t,obj.p(:,1),obj.p(:,2),q);
                    else
                        trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                        %trisurf(obj.t,obj.p(:,1),obj.p(:,2),q)
                        view(2); %shading interp;
                    end
                    cmocean('deep',numticks-1); cb = colorbar;
                    if logaxis
                        desiredTicks = round(10.^(linspace(min(q),max(q),numticks)),1);
                        caxis([log10(min(desiredTicks)) log10(max(desiredTicks))]);
                        cb.Ticks     = log10(desiredTicks);
                        for i = 1 : length(desiredTicks)
                            cb.TickLabels{i} = num2str(desiredTicks(i));
                        end
                    end
                    ylabel(cb,'m below geoid');
                    title('mesh topo-bathy');
                case('slp')
                    if proj
                        figure, h=m_trisurf(obj.t,obj.p(:,1),obj.p(:,2),...
                            hypot(obj.bx,obj.by)); view(2);
                    else
                        figure, h=trisurf(obj.t,obj.p(:,1),obj.p(:,2),...
                            hypot(obj.bx,obj.by),'facecolor', 'flat', 'edgecolor', 'none');
                        view(2);
                    end
                    colormap(cmocean('thermal'));
                    cb=colorbar; ylabel(cb,'slope');
                    caxis([0 0.25])
                case('ob') % outer boundary of mesh
                    [~,bpt] = extdom_edges2(obj.t,obj.p);
                    if proj
                       figure, m_plot(bpt(:,1),bpt(:,2),'r.');
                    else
                        figure, plot(bpt(:,1),bpt(:,2),'r.');
                    end
                case('reso')
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
                    if logaxis
                       q = log10(z); % plot on log scale with base
                    else
                       q = z; 
                    end
                    if proj
                        figure;
                        m_trimesh(obj.t,obj.p(:,1),obj.p(:,2),q);
                    else
                        figure;
                        trimesh(obj.t,obj.p(:,1),obj.p(:,2),q,'facecolor',...
                            'flat', 'edgecolor', 'none');
                        view(2); 
                    end
                    cmocean('thermal',numticks-1); cb = colorbar;
                    if logaxis
                        desiredTicks = round(10.^(linspace(min(q),max(q),10)));
                        caxis([log10(min(desiredTicks)) log10(max(desiredTicks))]);
                        cb.Ticks     = log10(desiredTicks);
                        for i = 1 : length(desiredTicks)
                            cb.TickLabels{i} = num2str(desiredTicks(i));
                        end
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
                        figure;
                        fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),values);
                        colormap(cmocean('deep',10));
                        caxis([0 5e-5])
                        colorbar;
                    else
                        display('Fort13 structure is empty!');
                    end
                case('cfvals')
                    if ~isempty(obj.f13)
                        ii = find(contains({obj.f13.defval.Atr(:).AttrName},'quadratic'));
                        defval  = obj.f13.defval.Atr(ii).Val;
                        userval = obj.f13.userval.Atr(ii).Val;
                        values = max(userval(2:end,:)',[],2);
                        figure;
                        fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),values);
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
                        figure;
                        fastscatter(obj.p(userval(1,:),1),obj.p(userval(1,:),2),values);
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
            if proj
                % now add the box
                m_grid('box','fancy','FontSize',12);
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
                disp('Renumbering the no flux boundary nodestrings...');
                for ib = 1 : obj.bd.nbou
                    tempcell{ib} = obj.bd.nbvv(1:obj.bd.nvell(ib),ib);
                end
                temp = cell2mat(tempcell');
                temp = perm_inv(temp)';
                temp = mat2cell(temp,cellfun(@length,tempcell));
                nbvv = zeros(size(obj.bd.nbvv,1),size(obj.bd.nbvv,2));
                for ib = 1 : obj.bd.nbou
                    for iv = 1 : obj.bd.nvell(ib)
                        nbvv(iv,ib) = temp{ib}(iv,:);
                    end
                end
                obj.bd.nbvv = nbvv;
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
        
        % make nodestrings
        function obj = makens(obj,type,dir)
            if nargin < 2
               error('Needs type: one of auto, islands or outer') 
            end
            trim = 0; periodic = 0;
            if strcmp(type(max(1,end-3):end),'trim')
                type = type(1:end-4); trim = 1;
            end
            switch type
                case('auto')
                    if ~isa(dir,'geodata')
                        error('third input must be a geodata class for auto makens')
                    else
                        gdat = dir;
                    end
                    
                    % Check for global mesh
                    L = find(obj.p(:,1) < -179, 1);
                    R = find(obj.p(:,1) > 179, 1);
                    if ~isempty(L) && ~isempty(R)
                        periodic = 1;
                        disp(['Detected global mesh applying automatic' ...
                              ' periodic BC fix'])
                    end
                    
                    % Get the boundaries
                    [etbv,~]  = extdom_edges2(obj.t,obj.p);
                    [~,poly_idx] = extdom_polygon(etbv,obj.p,1);
                    
                    % Get geodata outer and mainland polygons
                    outer = gdat.outer(~isnan(gdat.outer(:,1)),:);
                    nope = 0; neta = 0; nbou  = 0; nvel  = 0;
                    % Find the polygon that is a combination of ocean
                    % and mainland boundaries.
                    [~,idx] = sort(cellfun(@length,poly_idx),'descend');    
                    for op_ind = idx
                        % Get distance to outer
                        idv = poly_idx{op_ind};
                        [~, d_out] = ourKNNsearch(outer',obj.p(idv,:)',1);
                        if min(d_out) < gdat.h0/111e3
                            break;
                        end
                    end
                    % Get geodata outer and mainland polygons
                    outer = [gdat.outer(~isnan(gdat.outer(:,1)),:);
                             gdat.inner(~isnan(gdat.inner(:,1)),:)];
                    main = [gdat.mainland(~isnan(gdat.mainland(:,1)),:);
                            gdat.inner(~isnan(gdat.inner(:,1)),:)];
                    [~, d_out] = ourKNNsearch(outer',obj.p(idv,:)',1);    
                    [~, d_main] = ourKNNsearch(main',obj.p(idv,:)',1);
                    
                    % Mainland are nodes where shortest distance to
                    % mainland and outer are the same and where absolute
                    % distance to mainland is relatively small
                    if trim
                        mainland = abs(d_out - d_main) < gdat.h0/111e3 & ...
                                   d_main < 5*gdat.h0/111e3;
                    else
                        mainland = abs(d_out - d_main) < gdat.h0/111e3; 
                    end
                    
                    % indices of switch
                    Cuts  = find(diff(mainland ~= 0));
                    
                    if ~periodic
                        % Do not include open boundary that is smaller than 
                        % 10 vertices across
                        Cuts(diff(Cuts) < 10) = [];
                    end
                    
                    % Get length of largest island
                    % Delete the open boundary/mainland polygon
                    poly_idx(op_ind)= [];
                    L = max(1e3,max(cellfun(@length,poly_idx)));  
                    
                    if isempty(Cuts)
                        % if no mainland..
                        % Get ocean
                        nope = nope + 1;
                        nvdll(nope) = length(idv);
                        neta = neta + nvdll(nope);
                        ibtypee(nope) = 0;
                        nbdv(1:nvdll(nope),nope) = idv;
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
                                if ~mainland(Cuts(loop)+1)
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
                    
                    % loop through the remaining polygons (islands)
                    for poly_count = 1 : length(poly_idx)
                        idv = poly_idx{poly_count};
                        % islands
                        nbou = nbou + 1;
                        nvell(nbou) = length(idv);
                        nvel = nvel + nvell(nbou);
                        nbvv(1:nvell(nbou),nbou) = idv';
                        ibtype(nbou) = 21;
                    end
                    
                    % ocean boundary
                    obj.op.nope = nope ;
                    obj.op.neta = neta ;
                    obj.op.nvdll = nvdll ;
                    obj.op.ibtype = ibtypee ;
                    obj.op.nbdv = nbdv;
                    
                    % land boundary
                    obj.bd.nbou = nbou ;
                    obj.bd.nvel = nvel ;
                    obj.bd.nvell = nvell ;
                    obj.bd.ibtype = ibtype ;
                    obj.bd.nbvv = nbvv ;
                    
                    if ~periodic; return; end
                    
                    %% For periodic Bcs below
                    % Get open boundary points
                    opp = obj.op.nbdv(:); opp(opp == 0) = []; 
                    % Get the open boundary points on left side
                    I = obj.p(opp,1) < -179;
                    % ensure these are set to -180 correctly
                    obj.p(opp(I),1) = -180;
                    % use these as fixed points
                    % on left side...
                    pfixL = obj.p(opp(I),:);
                    eeL = Get_poly_edges([pfixL; NaN NaN]);
                    eeL(end,:) = [];
                    % apply flipped to right side..
                    pfixR = pfixL.* [-1, 1];
                    eeR = Get_poly_edges([pfixR; NaN NaN]);
                    eeR(end,:) = [];
                    pfix = [pfixL; pfixR];
                    ee = [eeL; eeR+length(pfixL)];
                    % Make dummy meshopts with the geodata
                    mshopts = meshgen();
                    mshopts.pfix = pfix;
                    mshopts.egfix = ee;

                    % Make a new triangulation using obj.p;
                    DT = delaunayTriangulation(obj.p);
                    % Now delete off right hand side and add on left hand side flipped
                    I  = find(DT.Points(:,1) > 179); DT.Points(I,:) = [];
                    DT.Points = [DT.Points; ...
                           DT.Points(DT.Points(:,1) < -179,:) .* [-1, 1]];

                    % Make polygon to delete internal triangles
                    bnde=extdom_edges2(obj.t,obj.p);
                    poly=extdom_polygon(bnde,obj.p,-1);
                    poly=cell2mat(poly');
                    m2 = obj; m2.p = DT.Points; m2.t = DT.ConnectivityList;
                    I1 = ourKNNsearch(m2.p(m2.t(:,1),:)',pfix',1);
                    I2 = ourKNNsearch(m2.p(m2.t(:,2),:)',pfix',1);
                    I3 = ourKNNsearch(m2.p(m2.t(:,3),:)',pfix',1);
                    bc = baryc(m2);
                    ee = Get_poly_edges(poly);
                    in = inpoly(bc,poly,ee);
                    % make sure fpix triangles remain
                    in(I1) = 1; in(I2) = 1; in(I3) = 1;
                    m2.t(~in,:) = [];
                    m2.b = []; m2.op = []; m2.bd = [];
                    mshopts.grd = m2; % get out the msh object
                    mshopts = clean(mshopts);

                    % Get the new match and cleaned grid
                    obj = mshopts.grd;

                    % Now add in the periodic point list
                    [I,d] = ourKNNsearch(obj.p',obj.p'+[360;0],1);
                    J = find(d < 1e-5);
                    % list of points that are same on both sides
                    periodic_bc_list = [I(J) J];    
                    obj.bd.nbou = 1 ;
                    obj.bd.nvel = length(periodic_bc_list) ;
                    obj.bd.nvell = obj.bd.nvel ;
                    obj.bd.ibtype = 94 ;
                    obj.bd.nbvv = periodic_bc_list ;

                case('islands')
                    [etbv,~]  = extdom_edges2(obj.t,obj.p);
                    [poly,poly_idx,max_ind] = extdom_polygon(etbv,obj.p,1);
                    
                    nbou = 0;
                    nvel = 0;
                    % the largest polygon will be a combination of ocean and mainland
                    % boundaries. Deal with this first, then remove it from the polygon
                    poly(max_ind) = [];
                    poly_idx(max_ind)=[];
                    
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
                    obj.bd.nbou = nbou ;
                    obj.bd.nvel = nvel ;
                    obj.bd.nvell = nvell ;
                    obj.bd.ibtype = ibtype ;
                    obj.bd.nbvv = nbvv ;
                    
                case('outer')
                    [bnde,bpts]=extdom_edges2(obj.t,obj.p);
                    
                    % use this to figure out the vstart and vend
                    figure, plot(bpts(:,1),bpts(:,2),'k.');
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
                    % Only add the periodic point list
                    [I,d] = ourKNNsearch(obj.p',obj.p'+[360;0],1);
                    J = find(d < 1e-5);
                    % list of points that are same on both sides
                    periodic_bc_list = [I(J) J];    
                    obj.bd.nbou = 1 ;
                    obj.bd.nvel = length(periodic_bc_list) ;
                    obj.bd.nvell = obj.bd.nvel ;
                    obj.bd.ibtype = 94 ;
                    obj.bd.nbvv = periodic_bc_list ;
                return;
                    
            end
            function txt = myupdatefcn2(~,event_obj,myarray)
                pos = get(event_obj,'Position');
                ind = find(abs(myarray(:,1)-pos(1))<eps & abs(myarray(:,2)-pos(2))<eps);
                txt = {['X: ',num2str(pos(1))],...
                    ['Y: ',num2str(pos(2))],...
                    ['Index: ',num2str(ind')]};
            end
        end
        
        function obj3=plus(obj1,obj2,varargin)
            % Merge 2-d simplex meshes compromised of np vertices and nt triangles
            % contained in the msh obejcts obj1 and obj2. Uses matlab's
            % implementation of the Boywer-Watson incremental triangulation.
            
            % The meshes must overlap partially (0.25 degrees is
            % recommended)
            
            % INPUTS:
            % mesh1: msh() class of base mesh.
            % mesh2: msh() class of mesh to be merged into the base mesh.
            
            % OUPUTS:
            % a msh object with...
            % pm: np x 2 coordinates of vertices of the merged mesh.
            % tm: nt x 3 matrix representing the 2-d simplices.
            % kjr, und, chl, sept. 2017 Version 1.0.
            if nargin > 2
                % for domains that are from ExtractSubdomain
                scheme = 1;
            else
                scheme = 0;
            end
            p1=obj1.p;
            t1=obj1.t;
            % form outer polygon 1
            disp('Forming outer boundary for base mesh...')
            bnde=extdom_edges2(t1,p1);
            poly1=extdom_polygon(bnde,p1,1);
            for i = 1 : length(poly1)
                poly1{i} = [poly1{i} ; NaN NaN];
            end
            poly_vec1=cell2mat(poly1');
            
            p2=obj2.p;
            t2=obj2.t;
            % form outer polygon 2
            disp('Forming outer boundary for mesh to be merged into base mesh...')
            bnde=extdom_edges2(t2,p2);
            poly2=extdom_polygon(bnde,p2,1);
            for i = 1 : length(poly2)
                poly2{i} = [poly2{i} ; NaN NaN];
            end
            poly_vec2=cell2mat(poly2');
            
            % delete the region in mesh 2 that is in the intersection with mesh 1.
            disp('Calculating the intersection between the meshes...');
            [xa, ya] = polybool('intersection', poly_vec1(:,1), poly_vec1(:,2), poly_vec2(:,1), poly_vec2(:,2));
            [edges]=Get_poly_edges([xa,ya; NaN NaN]);
            in1=inpoly(p2(t2(:,1),:),[xa,ya],edges);
            in2=inpoly(p2(t2(:,2),:),[xa,ya],edges);
            in3=inpoly(p2(t2(:,3),:),[xa,ya],edges);
            t2(in1 | in2 | in3,:)=[];
            p2=unique(p2(t2(:),:),'rows');
            
            
            disp('Cleaning up the intersection region between the meshes...');
            % figure out which triangles in t1 are on the boundary and
            % compute the circumcenters for them.
            [VToE]=VertToEle(t1);
            [~,on]=inpoly(p1,[xa,ya],edges);
            connT = VToE(:,on);
            connT = connT(:);
            connT(connT==0,:)=[];
            % for these boundary triangles compute the circumcenters.
            % compute the tighest fitting alpha shape around these
            % circumcenters. Compute the boundary of this alpha shape and
            % remove any points in p2 that are in the alphashape.
            TR1 = triangulation(t1,p1(:,1),p1(:,2));
            [cc,rc]=circumcenter(TR1,connT);
            % make sure all cc are in the neighborhood of
            % the intersection area.
            minx = min(xa); maxx = max(xa);
            miny = min(ya); maxy = max(ya);
            iboubox = [minx maxy
                minx miny
                maxx miny
                maxx maxy
                minx maxy];
            iboubox(:,1) = 2*iboubox(:,1)+(1-2)*nanmean(iboubox(1:end-1,1));
            iboubox(:,2) = 2*iboubox(:,2)+(1-2)*nanmean(iboubox(1:end-1,2));
            in = inpoly(cc,iboubox);
            % delete circumcenters that are supurious
            cc(~in,:) = []; rc(~in,:) = [];
            [tempy,tempx] = scircle1(cc(:,2),cc(:,1),rc);
            circpts = [tempx(:),tempy(:)];
            shp = alphaShape(circpts(:,1),circpts(:,2));
            kk = convhull(circpts(:,1),circpts(:,2));
            shppts = shp.Points;
            k = boundary(shppts(:,1),shppts(:,2));
            % delete p2 points that are in the cc's of the bndry
            % t1 triangles.
            in =inpoly(p2,[shppts(k,1),shppts(k,2)]);
            p2(in,:) = [];
            
            % now merge two meshes p1 and p2 to incrementally modify the
            % triangulation.
            disp('Merging meshes...')
            DTbase = delaunayTriangulation(p1(:,1),p1(:,2));
            DTbase.Points(end+(1:length(p2)),:)=p2;
            
            % now we won't retriangulate again, so unload the data
            pm=DTbase.Points; tm=DTbase.ConnectivityList;
            
            % prune triangles outside both domains.
            disp('Pruning outer triangles...')
            pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/3;
            
            % in1 is inside the first boundary polygon
            [edges1]=Get_poly_edges(poly_vec1);
            in1=inpoly(pmid,poly_vec1,edges1);
            
            % in2 is inside the second boundary polygon
            [edges2]=Get_poly_edges(poly_vec2);
            in2=inpoly(pmid,poly_vec2,edges2);
            
            minx = min(xa); maxx = max(xa);
            miny = min(ya); maxy = max(ya);
            iboubox = [minx maxy
                minx miny
                maxx miny
                maxx maxy
                minx maxy];
            % in3 is in intersection region. intersection region may not be
            % rectnagular.
            
            in3=inpoly(pmid,[iboubox(:,1),iboubox(:,2)]);
            
            % keep geographical "set difference" of the two meshes.
            if(scheme==1)
                tm((~in1 & ~in2),:)=[];
            else
                tm((~in1 & ~in2) | (in3 & ~in1),:)=[];
            end
            
            % clean up some more to avoid non-unique boundary edges.
            %disp('Fixing boundary problems...');
            [pm,tm]=Make_Mesh_Boundaries_Traversable(pm,tm,1,0.1);
            
            disp('Smoothing intersection...')
            % kjr 201803, smooth only the interesection region of the
            % meshes using an implicit smoother.
            pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/3;
            xt = circpts(kk,1);
            yt = circpts(kk,2);
            in = inpoly(pmid,[xt,yt]);
            tsub = tm(in,:); nnode = unique(tsub(:));
            psub = pm(nnode,:);
            % local-to-global mapping
            for i = 1 : length(psub)
                l2g(i) = nnode(i);
            end
            % global-to-local mapping
            for i = 1 : length(psub)
                g2l(nnode(i)) = i;
            end
            % renumber tsub
            for i = 1 : length(tsub)
                for j = 1 : 3
                    tsub(i,j) = g2l(tsub(i,j));
                end
            end
            [psub,tsub]=smoothmesh(psub,tsub,[],150,0.01);
            [psub,tsub]=direct_smoother_lur(psub,tsub,[],1);
            for i = 1 : length(psub)
                pm(l2g(i),:) = psub(i,:);
            end
            
            obj3=msh();
            obj3.p=pm; obj3.t=tm;
            
            %figure; simpplot(pm,tm);
            %hold on; plot(xa,ya,'r--')
            %hold on, plot(poly_vec1(:,1),poly_vec1(:,2),'k');
            %hold on; plot(bpts1(:,1),bpts1(:,2),'r.');
            %title('merged');
            
        end
        
        function out1 = CalcCFL(obj,dt)
            g      = 9.81;        % gravity
            R      = 6.3782064e6; % mean radius of earth
            % Convert lat-lon to radians
            pp = deg2rad(obj.p);
            % Get centre point of CPP projection
            Lon_C = mean(pp(:,1)); Lat_C = mean(pp(:,2));
            % Convert latitude-longtitude to metres
            X(:,1) = R*(pp(:,1) - Lon_C)*cos(Lat_C);
            X(:,2) = R*pp(:,2);
            % Get nearest two neighbours
            [idx,~] = ourKNNsearch(X',X',2) ;
            idx2 = idx(:)*0;
            % Rearrange to vector
            for ii = 1:length(idx)
               idx2(2*ii-1:2*ii) = idx(ii,:)'; 
            end
            d = m_lldist(obj.p(idx2,1),obj.p(idx2,2));
            d = d(1:2:end)*1e3;
            
            %U = sqrt(g*max(obj.b,1));
            %U(obj.b <= 0) = 1; % <--assume 1 m/s overland velocity.
            
            % wavespeed in ocean (second term represents orbital
            % velocity at 0 degree phase for 1-m amp. wave).
            U = sqrt(g*max(obj.b,1)) + sqrt(g./max(obj.b,1));
            if nargin > 1
                % Get CFL from input dt
                CFL = dt*U./d;  % <-- from the wave celerity
                out1 = CFL;
            else
                CFL = 1.0;
                dt = d.*CFL./U; % <-- from the wave celerity
                out1 = dt;
            end
        end
        
        
        function obj = CheckTimestep(obj,dt,varargin)
            %% Decimate mesh to achieve CFL condition for stability.
            % Takes a mesh and removes triangles and nodes to produce a mesh
            % that satisfies the given timestep requirements of the user by
            % incrementally modifying the triangulation nearby each edge
            % that violates the CFL.
            % kjr, chl, und, 2017
            %%
            desCFL     = 0.50;    %<-- set desired cfl (generally less than 0.80 for stable).
            if nargin < 3
                desIt = inf;        %<-- desired number of iterations
            else
                desIt = varargin{1};
            end
            %%
            F = scatteredInterpolant(obj.p(:,1),obj.p(:,2),obj.b,'linear','none');
            %%
            it  = 0;
            CFL = 999;
            tic
            while 1
                CFL = CalcCFL(obj,dt);
                bad = real(CFL) > desCFL;
                display(['Number of CFL violations ',num2str(sum(bad))]);
                disp(['Max CFL is : ',num2str(max(real(CFL)))]);
                if it==desIt,   break; end
                if sum(bad)==0 , break; end
                obj = DecimateTria(obj,bad);
                obj.b = F(obj.p(:,1),obj.p(:,2));
                it = it + 1;
            end
            toc
            disp(['Achieved max CFL of ',num2str(max(real(CFL))),' after ',num2str(it),' iterations.']);
            disp('Remove poor quality elements and fix connecitivity problems..');
            
            % Just delete really poor elements
            %[obj.p,obj.t] = fixmesh(obj.p,obj.t);
            %tq = gettrimeshquan(obj.p,obj.t);
            %obj.t(abs(tq.qm) < 0.10,:) = [];
            obj = Make_Mesh_Boundaries_Traversable( obj, 0.25, 0 );
            % Ensuring good numbering
            [obj.p,obj.t] = fixmesh(obj.p,obj.t);
            
            % add bathy back on
            obj.b = F(obj.p(:,1),obj.p(:,2));
            
            disp('Deleting boundary information...please renumber mesh');
            obj.bd = []; obj.op = [];
            return;
            
            function obj = DecimateTria(obj,bad)
                obj = Make_Mesh_Boundaries_Traversable(obj,0.001,1);
                % form outer polygon of mesh for cleaning up.
                bnde=extdom_edges2(obj.t,obj.p);
                poly1=extdom_polygon(bnde,obj.p,1);
                for i = 1 : length(poly1)
                    poly1{i} = [poly1{i} ; NaN NaN];
                end
                poly_vec1=cell2mat(poly1');
                [edges1]=Get_poly_edges(poly_vec1);
                % delete the points that violate the cfl, locally retriangulate, and then locally smooth
                DTbase = delaunayTriangulation(obj.p(:,1),obj.p(:,2));
                DTbase.Points(find(bad),:) = [];
                pm = DTbase.Points; tm = DTbase.ConnectivityList;
                pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/3;
                in1=inpoly(pmid,poly_vec1,edges1);
                tm(~in1,:)=[];
                [pm,tm]=fixmesh(pm,tm);
                if sum(bad) < 100e3
                    % find all the points nearby each "bad" point.
                    idx = ourKNNsearch(pm',obj.p(find(bad),:)',12) ;
                    idx = idx(:);
                    idx = unique(idx);
                    constr=setdiff((1:length(pm))',idx);
                    [pm,tm]=smoothmesh(pm,tm,constr,50,0.01);
                else
                    error('Adapation would result in potentially catastrophic lose of connectivity, try adapting to a smaller timestep');
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
            obj.p = []; obj.b = []; obj.t =[];
            obj.p = [fem_struct.x,fem_struct.y];
            obj.b = fem_struct.z;
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
        
        
        
    end % end methods
    
    
    
end % end class
