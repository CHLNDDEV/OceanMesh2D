classdef edgefx
    %   EDGEFX: Edgefunction class
    %   Constructs edgefunctions that are based on numerous geomteric and
    %   topographic criteria that guide the spatially variability of elemental
    %   resolution when building a mesh.
    %   Copyright (C) 2018  Keith Roberts & William Pringle
    %
    %   The following properties (with default values) are available for input
    %   to the edgefx method:
    %      defval = 0; % placeholder value if arg is not passed.
    %      addOptional(p,'dis',defval);
    %      addOptional(p,'fs',defval);
    %      addOptional(p,'wl',defval);
    %      addOptional(p,'slp',defval);
    %      addOptional(p,'ch',defval);
    %      addOptional(p,'min_el_ch',100);
    %      addOptional(p,'AngOfRe',60);
    %      addOptional(p,'max_el',inf);
    %      addOptional(p,'max_el_ns',inf);
    %      addOptional(p,'g',0.20);
    %      addOptional(p,'geodata',defval)
    %      addOptional(p,'lmsl',defval)
    %      addOptional(p,'dt',-1);
    %      addOptional(p,'fl',defval);
    %      addOptional(p,'Channels',defval);
    %      addOptional(p,'h0',defval);
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
        nx % number of x coords
        ny % number of y coords
        x0y0 % 2-tuple of bot. left of domain

        lmsl % lmsl boundary as a geodata class

        bbox
        boubox
        used % edge function keywords

        dis  % percent the edgelength changes in space
        fd   % distance function
        hhd  % matrix containing distance fx values.
        fs   % number of elements to resolve a feature
        fsd  % matrix containing feature size fx values.
        wl   % wavelength
        wld  % matrix containing wavelength fx values.
        slp  % slope edge function
        slpd % matrix containing filtered slope fx values.
        fl   % slope filter parameters
        ch   % channels
        chd  % matrix containing channel fx values.
        AngOfRe % Angle of reslope
        g    % max allowable grade of mesh.
        dt   % theoretical simulateable timestep

        F    % edge function in gridded interpolant format.

        h0    % min. resolution of edgelength function in meters
        gridspace % edgelength function resolution in WGS84 degrees.
        max_el % max. resolution in domain
        max_el_ns % max. resolution +-<0.01 from shoreline
        boudist % distance in WGS84 degrees to boundary of mesh

        min_el_ch % min. element size along the channel
        Channels % cell array of streams.

    end

    methods

        %% Class constructor/pass it a shape class and edgefx params.
        function obj = edgefx(varargin)
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

            p = inputParser;

            defval = 0; % placeholder value if arg is not passed.
            % add name/value pairs
            addOptional(p,'dis',defval);
            addOptional(p,'fs',defval);
            addOptional(p,'wl',defval);
            addOptional(p,'slp',defval);
            addOptional(p,'ch',defval);
            addOptional(p,'min_el_ch',100);
            addOptional(p,'AngOfRe',60);
            addOptional(p,'max_el',inf);
            addOptional(p,'max_el_ns',inf);
            addOptional(p,'g',0.20);
            addOptional(p,'geodata',defval)
            addOptional(p,'lmsl',defval)
            addOptional(p,'dt',-1);
            addOptional(p,'fl',defval);
            addOptional(p,'Channels',defval);
            addOptional(p,'h0',defval);

            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp = p.Results;
            % get the fieldnames of the edge functions
            inp = orderfields(inp,{'max_el','min_el_ch','AngOfRe',...
                'geodata','lmsl','Channels',...
                'dis','fs','fl','g','max_el_ns',...
                'wl','slp','ch','dt','h0'});
            fields = fieldnames(inp);
            % loop through and determine which args were passed.
            % also, assign reasonable default values if some options were

            % not assigned.
            for i = 1 : numel(fields)
                type = fields{i};
                switch type
                    % parse aux options first
                    case('max_el')
                        obj.max_el = inp.(fields{i});
                        if ~isempty(obj.max_el)
                            obj.max_el = inp.(fields{i});
                        end
                    case('max_el_ns')
                        obj.max_el_ns= inp.(fields{i});
                        if obj.max_el_ns ~=0
                            obj.max_el_ns = inp.(fields{i});
                        end
                        % obj.max_el_ns = obj.max_el_ns/111e3;
                    case('g')
                        obj.g= inp.(fields{i});
                        if obj.g ~=0
                            obj.g = inp.(fields{i});
                        end
                    case('fl')
                        obj.fl = inp.(fields{i});
                        if obj.fl ~=0
                            obj.fl = inp.(fields{i});
                        end
                    case('geodata')
                        if isa(inp.(fields{i}),'geodata')
                            feat = inp.(fields{i});
                        end
                    case('lmsl')
                        if isa(inp.(fields{i}),'geodata')
                            obj.lmsl = inp.(fields{i});
                        end
                    case('dt')
                        obj.dt=inp.(fields{i});
                        if obj.dt >= 0
                            obj.dt = inp.(fields{i});
                        end
                    case('h0')
                        obj.h0=inp.(fields{i});
                        if obj.h0 ~= 0
                            obj.h0 = inp.(fields{i});
                        else
                            obj.h0 = feat.h0;
                        end
                    case('min_el_ch')
                        obj.min_el_ch=inp.(fields{i});
                        if obj.min_el_ch ~= 100
                            obj.min_el_ch = inp.(fields{i});
                        end
                    case('AngOfRe')
                        obj.AngOfRe=inp.(fields{i});
                        if obj.AngOfRe ==60
                            obj.AngOfRe = inp.(fields{i});
                        end
                    case('Channels')
                        obj.Channels=inp.(fields{i});
                        if ~isempty(obj.Channels) ~=0
                            obj.Channels = inp.(fields{i});
                        end
                end
            end

            % kjr april 28, 2018-form mesh size grid on-the-fly
            obj.fd       = @dpoly;
            obj.x0y0     = feat.x0y0 + sqrt(eps);
            obj.gridspace = obj.h0/111e3;
            obj.nx       = ceil(abs(feat.x0y0(1)-feat.bbox(1,2))/obj.gridspace);
            obj.ny       = ceil(abs(feat.x0y0(2)-feat.bbox(2,2))/obj.gridspace);
            obj.bbox     = feat.bbox;
            obj.boubox   = feat.boubox;

            % WJP: Do Mercator projection for calculating
            % distances from shoreline
            if feat.bbox(1,2) > 180
                m_proj('mercator','lat', [-88 90], 'long',[0 360]);
            else
                m_proj('mercator','lat', [-88 90],'long',[-180 180]);
            end
            % now turn on the edge functions
            for i = 1 : numel(fields)
                type = fields{i};
                switch type
                    case('dis')
                        obj.dis = inp.(fields{i});
                        if obj.dis ~= 0
                            disp('Building distance function...');
                            obj = distfx(obj,feat);
                            obj.used{end+1} = 'dis';
                        end
                    case('fs')
                        obj.fs  = inp.(fields{i});
                        if obj.fs ~= 0
                            disp('Building feature size function...');
                            try
                                obj = featfx(obj,feat);
                                obj.used{end+1} = 'fs';
                            catch
                                disp('man, no medial points!')
                            end
                        end
                    case('wl')
                        obj.wl  = inp.(fields{i});
                        if obj.wl(1)~=0 && ~isempty(feat.Fb)
                            disp('Building wavelength function...');
                            obj = wlfx(obj,feat);
                            obj.used{end+1} = 'wl';
                        end
                    case('slp')
                        obj.slp  = inp.(fields{i});
                        if obj.slp(1) ~= 0 && ~isempty(feat.Fb)
                            disp('Building slope function...');
                            obj = slpfx(obj,feat);
                            obj.used{end+1} = 'slp';
                        end
                    case('ch')
                        obj.ch  = inp.(fields{i});
                        if obj.ch ~= 0
                            disp('Building channel function...');
                            obj = chfx(obj,feat);
                            obj.used{end+1} = 'ch';
                        end
                    case{'g','geodata','lmsl','max_el','min_el_ch',...
                            'Channels','max_el_ns','h0','dt','fl','AngOfRe'}
                        % dummy to avoid warning
                    otherwise
                        warning('Unexpected edge function name/value pairs.')
                end
                if i==numel(fields)
                    obj = finalize(obj,feat);
                    disp('Finalized edge function!');
                end
            end
        end

        %% Traditional distance function, linear distance from coastline.
        function obj = distfx(obj,feat)
            [d,obj] = get_dis(obj,feat);
            obj.boudist = d;
            obj.hhd = obj.h0 + obj.dis*abs(d);
        end
        %% Feature size function, approximates width of nearshore geo.
        function obj = featfx(obj,feat)
            if isempty(feat.mainland) && isempty(feat.inner)
                % No medial points so use distance function using grade
                warning('No mainland or inner')
                xg = CreateStructGrid(obj);
                obj.fsd = xg*NaN;
                clearvars xg
                return
            end
            % Make sure we don't create a singularity on the coast in the
            % distance function!
            [d,obj] = get_dis(obj,feat);
            obj.boudist = d ;

            [xg,yg] = CreateStructGrid(obj);
            % use a harvestine assumption
            dx = obj.h0*cosd(yg(1,:)); % for gradient function
            dy = obj.h0;               % for gradient function
            % Calculate the gradient of the distance function.
            [ddy,ddx] = EarthGradient(d,dy,dx);
            d_fs = sqrt(ddx.^2 + ddy.^2);
            clearvars ddx ddy
            % WJP: This fix is to put a medial point in narrow channels
            rf = []; cf = [];
            for ii = 2:size(d,1)-1
                for jj = 2:size(d,2)-1
                    if d(ii,jj) < 0
                        if (d(ii-1,jj) >= 0 && ...
                                d(ii+1,jj) >= 0) || ...
                                (d(ii,jj-1) >= 0 && ...
                                d(ii,jj+1) >= 0)
                            rf(end+1) = ii;
                            cf(end+1) = jj;
                        end
                    end
                end
            end
            % Find singularties in the distance function that are
            % within the poly to get the medial axis.
            % Lets only create medial points in places that are
            % sufficiently far away from the coastline to capture the
            % feature and thus aren't spurious.
            d_fs = reshape(d_fs,[],1); d = reshape(d,[],1);
            lg = d_fs < 0.90  & d < -0.5*obj.h0 ;
            [r,c] = ind2sub([obj.nx obj.ny],find(lg));
            % create the coordinates using x0y0 + r*h0
            r = [r; rf']; c = [c; cf']; %WJP: add on "fixed points"
            x_kp = obj.x0y0(1) + (r-1)*obj.gridspace;
            y_kp = obj.x0y0(2) + (c-1)*obj.gridspace;
            clearvars lg r c cf rf
            % Ensure there are enough medial points nearby
            % to avoid surpurious medial points.
            % WJP: We want a line of points larger than around 7. This
            % corresponds to having three points back or forward up the
            % line. Let's check for those three closest points and ensure
            % they are within about co*h0 distance from each other.
            % co is the cutoff distance = 2*sqrt(2) (two diagonal dist)
            co = 2*sqrt(2);
            if length(x_kp) > 12
                [~, dmed] = WrapperForKsearch([x_kp,y_kp],[x_kp,y_kp],4);
                prune = dmed(:,2) > co*obj.h0 | dmed(:,3) > 2*co*obj.h0 | ...
                    dmed(:,4) > 3*co*obj.h0;
                x_kp( prune ) = [];
                y_kp( prune ) = [];
            end

            % Now get the feature size along the coastline
            if length(x_kp) > 12
                % Use KD-tree
                dPOS = 0*d;
                noblks = ceil(length(dPOS)*2*8*1e-9);
                blklen = floor(length(dPOS)/noblks);
                ns = 1;
                for blks = 1:noblks
                    if blks == noblks
                        ne = length(dPOS);
                    else
                        ne = ns + blklen - 1;
                    end
                    [~,dPOS(ns:ne)] = WrapperForKsearch([x_kp,y_kp],...
                        [xg(ns:ne)',yg(ns:ne)'],1);
                    ns = ne + 1;
                end
            else
                % No medial points so use distance function using grade
                warning('No medial points, resorting to distance function')
                d = reshape(d,obj.nx,[]);
                obj.fsd = obj.h0 + obj.g*abs(d);
                clear x_kp y_kp d d_fs dPOS xg yg
                return
            end
            % reshape back
            d = reshape(d,obj.nx,[]); dPOS = reshape(dPOS,obj.nx,[]);
            % Feature_size is distance from medial axis plus distance to
            % coastline. min_el is then feature_size*2/R where R is
            % number of elements to model the feature
            % WJP: d needs to be the absolute value because when we
            % include the floodplain since dPOS is not signed but d is,
            % they will cancel out producing really fine resolution
            W = dPOS+abs(d);
            if obj.fs < 0
                % automatic feature size
                fs_auto = min(-obj.fs,ceil(W/obj.h0));
                obj.fsd = 2*W./fs_auto;
            else
                obj.fsd = 2*W/obj.fs;
            end
            clear x_kp y_kp d d_fs dPOS xg yg
        end

        function [d,obj] = get_dis(obj,feat)
            % Function used by distfx and featfx to return distance
            % and make the Fdis interpolant
            d = feval(obj.fd,obj,feat);
            d = reshape(d,obj.nx,[]);
        end

        %% Wavelength edgefx.
        function obj = wlfx(obj,feat)

            % interpolate DEM's bathy linearly onto our edgefunction grid.
            [xg,yg] = CreateStructGrid(obj);
            tmpz    = feat.Fb(xg,yg);
            if ~isempty(feat.Fb2)
               tmpz2 = feat.Fb2(xg,yg);
               tmpz(isnan(tmpz)) = tmpz2(isnan(tmpz));
               clear tmpz2
            end
            % Initialise wld
            obj.wld = NaN([obj.nx,obj.ny]);
            grav    = 9.807;
            period  = 12.42*3600; % M2 period in seconds

            for param = obj.wl'
                if numel(param)==1
                    % no bounds specified.
                    wlp = param(1);
                    % set cuttof at 10 m by default
                    dp1 = -10;
                    dp2 = -inf;
                else
                    wlp = param(1);
                    dp1 = param(2);
                    dp2 = param(3);
                end
                % limit to 1 m
                twld = period*sqrt(grav*max(abs(tmpz),1))/wlp;
                % Set wld with mask applied
                obj.wld(tmpz < dp1 & tmpz > dp2 ) = ...
                    twld(tmpz < dp1 & tmpz > dp2);
                clearvars twld
            end
            clearvars tmpz xg yg;
        end

        %% Topographic length scale/slope edge function.
      %% Topographic length scale/slope edge function.
        function obj = slpfx(obj,feat)

            [xg,yg] = CreateStructGrid(obj);

            tmpz    = feat.Fb(xg,yg);
            if ~isempty(feat.Fb2)
               tmpz2 = feat.Fb2(xg,yg);
               tmpz(isnan(tmpz)) = tmpz2(isnan(tmpz));
               clear tmpz2
            end
            tmpz(tmpz > 50) = 50; % ensure no larger than 50 m above land
            % use a harvestine assumption
            dx = obj.h0*cosd(min(yg(1,:),85)); % for gradient function
            dy = obj.h0;               % for gradient function
            % lets filter the bathy to get only relevant features
            % loop over each set of bandpass filter lengths
            tmpz_f = zeros(size(tmpz));
            if obj.fl(1) < 0 && obj.fl(1) ~= -999
                disp('INFO: Rossby radius of deformation filter on.') ;
                rbfilt = abs(obj.fl(1));
                barot  = 1;
                if length(obj.fl) > 1 && obj.fl(2) < 0
                    disp('INFO: Using 1st-mode baroclinic Rossby radius.');
                    barot = 0;
                else
                    disp('INFO: Using barotropic Rossby radius.') ;
                end
                obj.fl = [];
                filtit = 1;
            elseif obj.fl(1) == 0
                disp('INFO: Slope filter is off.');
                obj.fl = [];
                tmpz_f = tmpz;
                filtit = 0 ;
            elseif obj.fl(1) == -999
                disp('INFO: Slope filter is LEGACY.');
                obj.fl = [];
                tmpz_f = tmpz;
                filtit = -999;
            else
                filtit = 0;
                for lambda = obj.fl'
                    if length(lambda) == 1 || lambda(2) == 0
                        % do a low pass filter
                        tmpz_ft  = filt2(tmpz,dy,lambda(1),'lp') ;
                    elseif all(lambda ~= 0)
                        % do a bandpass filter
                        tmpz_ft  = filt2(tmpz,dy,lambda,'bp') ;
                    else
                        % highpass filter not recommended
                        warning(['Highpass filter on bathymetry in slope' ...
                            'edgelength function is not recommended'])
                        tmpz_ft  = filt2(tmpz,dy,lambda(2),'hp') ;
                    end
                    tmpz_f = tmpz_f + tmpz_ft;
                end
            end

            % Rossby radius of deformation filter
            if filtit == 1
                tic
                bs = NaN([obj.nx,obj.ny]);
                % break into 10 deg latitude chunks, or less if higher res
                div = ceil(min(1e7/obj.nx,10*obj.ny/(max(yg(:))-min(yg(:)))));
                grav = 9.807;
                nb = ceil(obj.ny/div);
                n2s = 1;
                for jj = 1:nb
                    n2e = min(obj.ny,n2s + div - 1);
                    % Rossby radius of deformation filter
                    % See Shelton, D. B., et al. (1998): Geographical variability of the first-baroclinic Rossby radius of deformation. J. Phys. Oceanogr., 28, 433-460.
                    ygg = yg(:,n2s:n2e);
                    dxx = mean(dx(n2s:n2e));
                    f = 2*7.29e-5*abs(sind(ygg));
                    if barot
                        % barotropic case
                        c = sqrt(grav*max(1,-tmpz(:,n2s:n2e)));
                    else
                        % baroclinic case (estimate Nm to be 2.5e-3)
                        Nm = 2.5e-3;
                        c = Nm*max(1,-tmpz(:,n2s:n2e))/pi;
                    end
                    rosb = c./f;
                    clear f;
                    % update for equatorial regions
                    I = abs(ygg) < 5; Re = 6371e3;
                    twobeta = 4*7.29e-5*cosd(ygg(I))/Re;
                    rosb(I) = sqrt(c(I)./twobeta);
                    clear twobeta
                    % limit rossby radius to 10,000 km for practical purposes
                    rosb(rosb > 1e7) = 1e7;
                    % Keep lengthscales rbfilt * barotropic
                    % radius of deformation
                    rosb = min(10,max(0,floor(log2(rosb/dy/rbfilt))));
                    edges = double(unique(rosb(:)));
                    bst = rosb*0;
                    for i = 1:length(edges)
                        if edges(i) > 0
                            mult = 2^edges(i);
                            xl = 1; xu = obj.nx;
                            if (max(xg(:)) > 179 && min(xg(:)) < -179) || ...
                                    (max(xg(:)) > 359 && min(xg(:)) < 1)
                                % wraps around
                                xr = [obj.nx-mult/2+1:1:obj.nx xl:xu 1:mult/2];
                            else
                                xr = xl:xu;
                            end
                            yl = max(1,n2s-mult/2);
                            yu = min(obj.ny,n2e+mult/2);
                            if max(yg(:)) > 89 && yu == obj.ny
                                % create mirror around pole
                                yr = [yl:yu yu-1:-1:2*obj.ny-n2e-mult/2];
                            else
                                yr = yl:yu;
                            end
                            if mult == 2
                                tmpz_ft = filt2(tmpz(xr,yr),...
                                           [dxx dy],dy*2.01,'lp');
                            else
                                tmpz_ft = filt2(tmpz(xr,yr),...
                                           [dxx dy],dy*mult,'lp');
                            end
                            % delete the padded region
                            tmpz_ft(1:find(xr == 1,1,'first')-1,:) = [];
                            tmpz_ft(obj.nx+1:end,:) = [];
                            tmpz_ft(:,1:find(yr == n2s,1,'first')-1) = [];
                            tmpz_ft(:,n2e-n2s+2:end) = [];
                        else
                            tmpz_ft = tmpz(:,n2s:n2e);
                        end
                        [by,bx] = EarthGradient(tmpz_ft,dy,dx(n2s:n2e)); % get slope in x and y directions
                        tempbs  = sqrt(bx.^2 + by.^2);          % get overall slope
                        bst(rosb == edges(i)) = tempbs(rosb == edges(i));
                    end
                    bs(:,n2s:n2e) = bst;
                    n2s = n2e + 1;
                end
                toc
                % legacy filter
            elseif filtit==-999
                bs = NaN([obj.nx,obj.ny]);
                % Rossby radius of deformation filter
                f = 2*7.29e-5*abs(sind(yg));
                % limit to 1000 km
                rosb = min(1000e3,sqrt(9.81*abs(tmpz))./f);
                % autmatically divide into discrete bins
                [~,edges] = histcounts(rosb);
                tmpz_ft  = tmpz; dyb = dy;
                % get slope from filtered bathy for the segment only
                [by,bx] = EarthGradient(tmpz_ft,dy,dx); % get slope in x and y directions
                tempbs  = sqrt(bx.^2 + by.^2); % get overall slope
                for i = 1:length(edges)-1
                    sel = rosb >= edges(i) & rosb <= edges(i+1);
                    rosbylb = mean(edges(i:i+1));
                    if rosbylb > 2*dyb
                        disp(['i = ',num2str(i), ' rl/dx = ',num2str(rosbylb/dyb)])
                        tmpz_ft  = filt2(tmpz_ft,dyb,rosbylb,'lp');
                        dyb = rosbylb;
                        % get slope from filtered bathy for the segment only
                        [by,bx] = EarthGradient(tmpz_ft,dy,dx); % get slope in x and y directions
                        tempbs  = sqrt(bx.^2 + by.^2); % get overall slope
                    else
                        % otherwise just use the same tempbs from before
                    end
                    % put in the full one
                    bs(sel) = tempbs(sel);
                end
            else
                % get slope from (possibly filtered) bathy
                [by,bx] = EarthGradient(tmpz_f,dy,dx); % get slope in x and y directions
                bs      = sqrt(bx.^2 + by.^2); % get overall slope
            end
            clear bx by tmpz_f tmpz_ft

            % Allow user to specify depth ranges for slope parameter.
            obj.slpd = NaN([obj.nx,obj.ny]);
            for param = obj.slp'
                if numel(param)==1
                    % no bounds specified. valid in this range.
                    slpp = param(1);
                    % default cutoff is 10 m
                    dp1 = -10;
                    dp2 = -inf;
                else
                    slpp = param(1);
                    dp1 = param(2);
                    dp2 = param(3);
                end
                % Calculating the slope function
                dp = max(1,-tmpz);
                tslpd = (2*pi/slpp)*dp./(bs+eps);
                obj.slpd(tmpz < dp1 & tmpz > dp2 ) = ...
                    tslpd(tmpz < dp1 & tmpz > dp2);
                clearvars tslpd
            end
            clearvars tmpz xg yg
        end

        %% Channel edgefunction
        function obj = chfx(obj,feat)

            ang_of_reslope=obj.AngOfRe ;      % Estimate width of channel using tangent of this angle times depth.

            % STEP 1: Calculate the width of each channel using a v-shape approx to
            % channel's cross sectional area.
            for jj = 1:length(obj.Channels)
                if (isempty(obj.Channels{jj})); continue ;end
                pts{jj} = obj.Channels{jj};
                dp{jj} = feat.Fb(pts{jj});                               % depth at x,y channel locations
                radii{jj}=(tand(ang_of_reslope)*max(1,-dp{jj}));         % estimate of channel's width in degrees at x,y locations ASSUMING angle of reslope
                tempbb{jj} = feat.boubox;
            end

            [xg,yg] = CreateStructGrid(obj);

            % STEP 2: For each channel point, set the resolution around a
            % neighborhood of points as |h|/ch where ch is non-dim # between 0.5-1.5
            obj.chd = NaN([obj.nx,obj.ny]);
            % prune the potentially large channel network to only ones
            % partially enclosed in boubox
            isin = cellfun(@inpolysum,pts,tempbb);
            pts(isin==0) = []; radii(isin==0)= [];

            %figure; plot(xg,yg,'k.');
            tic
            tidx = [];
            for jj =1:length(pts)                                           % For each channel jj
                sel = pts{jj};
                in  = inpoly(sel,feat.boubox(1:end-1,:));
                sel(~in,:) = [];
                for jjj = 1 : size(sel,1)                                   % For each channel point jjj on channel jj.
                    % will the stencil be too large (nidx^2)?
                    nidx = ceil(radii{jj}(jjj)/obj.h0);                       % stencil of points around channel point.
                    % find linear index of channel point in 2-d matrix.
                    [lidx] = FindLinearIdx(sel(jjj,1),sel(jjj,2),xg,yg);
                    % convert to r,c or grid indices.
                    [r,c]     = ind2sub(size(xg),lidx);
                    [cid,rid] = ndgrid(c-nidx:c+nidx,r-nidx:r+nidx);    % grid of nearby points
                    % ensure that all these points are within the domain.
                    rid = min(max(rid(:),1),size(xg,1));
                    cid = min(max(cid(:),1),size(xg,2));
                    % convert back to linear indices.
                    temp = sub2ind(size(xg),rid,cid);
                    tidx = [tidx; temp];
                end
            end
            toc
            tmpz = feat.Fb(xg(tidx),yg(tidx));
            if ~isempty(feat.Fb2)
               tmpz2 = feat.Fb2(xg(tidx),yg(tidx));
               tmpz(isnan(tmpz)) = tmpz2(isnan(tmpz));
               clear tmpz2
            end
            dp = max(1,-tmpz);
            obj.chd(tidx) = dp/obj.ch;
            obj.chd(obj.chd < obj.min_el_ch) = obj.min_el_ch;
            clearvars Fb pts dp radii tempbb xg yg
        end


        %% General function to plot
        function plot(obj,type)
            % form grid here.
            % working out plottable size
            [xgrid,ygrid] = CreateStructGrid(obj);

            if nargin == 2
                switch type
                    case('dis')
                        val = obj.hhd;
                        tt = 'Distance EdgeLength Function';
                    case('fs')
                        val = obj.fsd;
                        tt = 'Feature Size EdgeLength Function';
                    case('wl')
                        val = obj.wld;
                        tt = 'Wavelength EdgeLength Function';
                    case('ch')
                        val = obj.chd;
                        tt = 'Channel EdgeLength Function';
                    case('slp')
                        val = obj.slpd;
                        tt = 'Slope EdgeLength Function';
                    otherwise
                        warning('Unexpected plot type. No plot created.')
                end
            else
                val = obj.F.Values;
                tt = 'Total EdgeLength Function';
            end

            % working out stride
            mem = inf; stride = obj.gridspace;
            while mem > 1
                xx = obj.x0y0(1):stride:obj.bbox(1,2);
                yy = obj.x0y0(2):stride:obj.bbox(2,2);
                xs = whos('xx'); ys = whos('yy');
                mem = xs.bytes*ys.bytes/1e9;
                stride = stride*2;
            end
            stride = ceil(stride/obj.gridspace);
            I = [1:stride:size(xgrid,1)];
            J = [1:stride:size(xgrid,2)];

            figure;
            m_proj('merc',...
                'long',[obj.bbox(1,1) obj.bbox(1,2)],...
                'lat',[obj.bbox(2,1) obj.bbox(2,2)])
            hold on; m_fastscatter(xgrid(I,J),ygrid(I,J),val(I,J));
            cb = colorbar; ylabel(cb,'edgelength in meters');
            colormap(lansey)
            m_grid() %'xtick',10,'tickdir','out',...
            % 'yaxislocation','left','fontsize',7);
            title(tt);
        end

        %% Finalize edge function
        % Add grading, cfl limiting and bound enforcement.
        function obj = finalize(obj,feat)
            % package the edge functions into a known order.
            counter = 0;
            hh = zeros([obj.nx,obj.ny]);
            for i = 1 : numel(obj.used)
                type = obj.used{i};
                switch type
                    case('dis')
                        counter = counter + 1;
                        hh(:,:,counter) = obj.hhd;
                        obj.hhd = single(obj.hhd);
                    case('fs')
                        counter = counter + 1;
                        hh(:,:,counter) = obj.fsd;
                        obj.fsd = single(obj.fsd);
                    case('wl')
                        counter = counter + 1;
                        hh(:,:,counter) = obj.wld;
                        obj.wld = single(obj.wld);
                    case('slp')
                        counter = counter + 1;
                        hh(:,:,counter) = obj.slpd;
                        obj.slpd = single(obj.slpd);
                    case('ch')
                        counter = counter + 1;
                        hh(:,:,counter) = obj.chd;
                        obj.ch = single(obj.chd);
                    otherwise
                        error('FATAL:  Could not finalize edge function');
                end
            end

            [hh_m] = min(hh,[],3);
            clearvars hh

            [xg,yg] = CreateStructGrid(obj);

            % Before grading, ensure edgefx reflects the spacing in
            % proximity to constraints.
            if ~isempty(feat.weirs)
                % form temporary mesh size function interpolat
                % if spacing is negative -999 then it will use the mesh
                % size function to determine the spacing.
                %Ftmp = griddedInterpolant(xg,yg,hh_m,'linear','nearest');

                for i =1:length(feat.weirs)
                    if iscell(feat.weirs)
                        % get spacing in meters along the face of the weir
                        weir_spacing = feat.weirs{i}(1,4);
                        % get crestline points in wgs84
                        txy = feat.weirs{i}(:,1:2) ;
                    else
                        % get spacing in meters along the face of the weir
                        weir_spacing = feat.weirs(i).min_ele;
                        % get crestline points in wgs84
                        txy = [feat.weirs(i).X feat.weirs(i).Y];
                    end
                    xy = [];
                    [xy(:,2),xy(:,1)] = my_interpm(txy(:,2),txy(:,1),...
                        (obj.h0/2)/111e3);
                    % edit the mesh size function in proximity to the weir
                    for j = 1 : length(xy)
                        nidx = 10 ;
                        % return the indices into the mesh size function
                        [lidx] = FindLinearIdx(xy(j,1),xy(j,2),xg,yg);
                        [r,c]  = ind2sub(size(xg),lidx);
                        [cid,rid]=ndgrid((c-nidx:c+nidx)',(r-nidx:r+nidx)');    % grid of nearby points
                        rid = [rid(:);r]; cid=[cid(:);c];
                        hh_m(rid,cid) = weir_spacing ;
                    end
                end
            end


            % enforce all mesh resolution bounds,grade and enforce the CFL in planar metres
            if ~isinf(obj.max_el_ns) && ~isempty(obj.boudist)
                nearshore = abs(obj.boudist) < 2/3*obj.h0;
                hh_m(nearshore & hh_m > obj.max_el_ns) = obj.max_el_ns;
            end
            hh_m(hh_m < obj.h0 ) = obj.h0;

            % Compute tmpz if necessary
            if length(obj.max_el) > 1 || length(obj.g) > 1 || obj.dt(1) >= 0
                tmpz    = feat.Fb(xg,yg);
                if ~isempty(feat.Fb2)
                   tmpz2 = feat.Fb2(xg,yg);
                   tmpz(isnan(tmpz)) = tmpz2(isnan(tmpz));
                   clear tmpz2
                end
            end

            for param = obj.max_el'
                if numel(param)==1 && param~=0
                    mx   = obj.max_el(1);
                    limidx = hh_m > mx | isnan(hh_m);

                    hh_m( limidx ) = mx;
                else
                    mx  = param(1);
                    dp1 = param(2);
                    dp2 = param(3);

                    limidx = (tmpz < dp1 & tmpz > dp2) & ...
                             (hh_m > mx | isnan(hh_m));

                    hh_m( limidx ) = mx;
                end
            end

            % Make sure this is called before releasing memory...
            if obj.dt(1) == 0
                % Find min allowable dt based on dis or fs function
                if any(~cellfun('isempty',strfind(obj.used,'dis')))
                    hh_d = obj.hhd;
                elseif any(~cellfun('isempty',strfind(obj.used,'fs')))
                    hh_d = obj.fsd;
                else
                    error(['FATAL: cannot use automatic timestep ' ...
                           'limiter without specifying dis or fs option']);
                end
            end

            obj = release_memory(obj) ;

            disp('Relaxing the mesh size gradient');
            if all(abs(obj.bbox(1,:)) == 180)
                % for global mesh make it cyclical
                hh_m = [hh_m(end,:); hh_m; hh_m(1,:)];
            end
            hfun = reshape(hh_m',[numel(hh_m),1]);

            dx = obj.h0*cosd(min(yg(1,:),85)); % for gradient function
            dy = obj.h0;                       % for gradient function

            % make g a function of space
            dmy     = xg*0 ;
            for param = obj.g'
                if numel(param)==1 && param~=0
                    lim   = obj.g(1);
                    dmy  = dmy + lim ;
                else
                    lim  = param(1);
                    dp1  = param(2);
                    dp2  = param(3);

                    limidx = (tmpz < dp1 & tmpz > dp2) ;

                    dmy( limidx ) = lim;
                end
            end
            if all(abs(obj.bbox(1,:)) == 180)
                % for global mesh make it cyclical
                dmy = [dmy(end,:); dmy; dmy(1,:)];
            end
            fdfdx = reshape(dmy',[numel(dmy),1]);
            clearvars dmy;

            [hfun,flag] = limgradStruct(obj.ny,dx,dy,hfun,...
                fdfdx,sqrt(length(hfun)));
            if flag == 1
                disp('Gradient relaxing converged!');
            else
                error(['FATAL: Gradient relaxing did not converge, '
                    'please check your edge functions']);
            end
            % reshape it back
            hh_m = reshape(hfun,obj.ny,[])';
            clearvars hfun fdfdx
            if all(abs(obj.bbox(1,:)) == 180)
                % for global mesh make it cyclical
                hh_m = hh_m(2:end-1,:);
                hh_m(end,:) = hh_m(1,:);
            end

            % Limit CFL if dt >= 0, dt = 0 finds dt automatically.
            if obj.dt(1) >= 0
                if isempty(feat.Fb)
                    error('No DEM supplied Can''t CFL limit.')
                end
                % Estimate wavespeed in ocean (second term represents orbital
                % velocity at 0 degree phase for 1-m amp. wave).
                grav = 9.807;
                % Limit the minimum depth to 1 m (ignore overland)
                tmpz(tmpz > - 1) = -1;
                u = sqrt(grav*abs(tmpz)) + sqrt(grav./abs(tmpz));

                % Desired timestep to bound edgefx
                desDt = obj.dt(1);
                if isnan(desDt)
                    disp('No timestep enforcement due to no distance function')
                end

                % Determine Courant min. and max. bounds
                if length(obj.dt) < 2
                    % default behavior if no bounds are passed
                    % for explicit type schemes we ensure the maxCr is
                    % bounded
                    minCr = eps; % Cannot be zero!!!
                    maxCr = 0.5;
                elseif length(obj.dt) == 2
                    if obj.dt(1)==0.0
                        error('dt = 0 is only intended for restricting the upper bound of Cr')
                    end
                    % minimum Courant number bound
                    minCr = obj.dt(2);
                    maxCr = inf;
                elseif length(obj.dt) == 3
                    % minimum and maximum Courant bound
                    minCr = obj.dt(2);
                    maxCr = obj.dt(3);
                else
                    error('Length of dt name-value pair should be <=3')
                end

                if minCr > maxCr
                    error('MinCr is > MaxCr, please switch order in dt name value call')
                end

                % Automatic timestep selection option (for ADCIRC models)
                if desDt == 0
                    hh_d(hh_d < obj.h0) = obj.h0;
                    obj.dt = min(min(maxCr*hh_d./u));
                    desDt  = obj.dt;
                end

                % Enforcement of Courant bounds in mesh size function
                fprintf(1, [ ...
                    'Enforcing timestep of ',num2str(desDt),' seconds ', ...
                    'with Courant number bounds of ',num2str(minCr),' to ',num2str(maxCr),'\n', ...
                    ] ) ;


                Cr = (desDt*u)./hh_m; % Courant number for given desDt for mesh size variations

                % resolution that meets max. Cr. criteria for desDt
                dxn_min = u*desDt/maxCr;
                % resolution that meets min. Cr. criteria for desDt
                dxn_max = u*desDt/minCr;

                % resolve min. Cr violations
                hh_m( Cr <= minCr) = dxn_max( Cr <= minCr); %--in planar metres

                % resolve max. Cr violations
                hh_m( Cr > maxCr) = dxn_min( Cr > maxCr);   %--in planar metres

                clear dxn_min dxn_max Cr
                clear u hh_d tmpz
            end


            obj.F = griddedInterpolant(xg,yg,hh_m,'linear','nearest');

            clearvars xg yg
        end

        function [xg,yg] = CreateStructGrid(obj)
            [xg,yg] = ndgrid(obj.x0y0(1) + (0:obj.nx-1)'*obj.gridspace, ...
                obj.x0y0(2) + (0:obj.ny-1)'*obj.gridspace);
        end

        function obj = release_memory(obj)
            % releases heavy components from data structures
            if isa(obj,'edgefx')
                disp('--------------------------------------------------');
                disp('Releasing individual edge functions from memory...');
                disp('--------------------------------------------------');
                if ~isempty(obj.fsd)
                    obj.fsd = [];
                end

                if ~isempty(obj.wld)
                    obj.wld = [];
                end

                if ~isempty(obj.slpd)
                    obj.slpd = [];
                end

                if ~isempty(obj.chd)
                    obj.chd = [];
                end

                if ~isempty(obj.ch)
                    obj.ch = [];
                end

                if ~isempty(obj.hhd)
                    obj.hhd = [];
                end

                if ~isempty(obj.boudist)
                    obj.boudist = [];
                end

                if ~isempty(obj.Channels)
                    obj.Channels = [];
                end

                if ~isempty(obj.ch)
                   obj.ch = [];
                end

                if ~isempty(obj.boudist)
                   obj.boudist = [];
                end
            end
        end


    end

end
