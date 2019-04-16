function obj = Make_f15( obj, ts, te, dt, varargin )
% obj = Make_f15(obj,ts, te, dt, varargin)
% Input a msh class object, the simulation start time, end times, timestep
% and any other optional arguments and get back a mesh class object with
% the f15 structure populated
% 
%  Inputs:   1) A msh class obj (requires nodestrings if you want to get
%               the elevation specified boundary forcings), requires f13
%               populated to put in the NWP names if you want any f13 options.
%            2) ts: the simulation start time in datestring format e.g.
%               ts = '01-Jan-2018 03:00'
%            3) te: the simulation end time in datestring format e.g.
%               te = '01-Feb-2018 22:00'
%            4) dt: the simulation time step in seconds
%            5) Optional name-value arguments:
%
%               'const' : the harmonic constituents to use for the tidal 
%               potential, elevation specified boundaries (if you have
%               nodestrings and a tidal database specified), and the 
%               harmonic analysis. Give value as a cell array of characters
%               e.g: C = {'M2','S2','K1','O1'}; or a string array,
%               e.g: C = ["M2","S2","K1","O1"]; or u can put:
%                    C = "major8" to do all the major eight constituents
%               Note that the tidal nodal factors and so forth will be
%               automatically calculated based on ts and te. 
%
%               'tidal_database' : the directory + filename tidal database of
%               elevations to interpolate to your boundaries. At the moment
%               just handles TPXO9.1 (Netcdf), available for download at:
%               http://volkov.oce.orst.edu/tides/global.html. The filename 
%               will be 'h_tpxo9.v1.nc'. This option requires the 'const' 
%               name-value pair input.
%
%               'sta_database' : a 1x2 cell, where the first element contains
%               a character referring to a specific database. At the moment
%               it only handles CO-OPS/NOS/NOAA stations. The second
%               element contains a vector of any of 1,2,3 corresponding to
%               elev, vel and met stations. e.g. set
%               ...('sta_database',{'CO-OPS',[1 2]}) to output the CO-OPS
%               stations for elevation and velocity recording. 
%
%  Outputs:  1) msh class obj with f15 struct populated
%
%  Author:      William Pringle and Keith Roberts                               
%  Created:     March 15 2018                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Putting in default options if f15 isn't already populated
if isempty(obj.f15) 
    % Run Description and ID
    f15dat.rundes = obj.title ;  % RUNDES
    f15dat.runid  = 'Run_001' ;  % RUNID
    % NFOVER
    f15dat.nfover = 1 ;
    % NABOUT
    f15dat.nabout = 0 ;
    % NSCREEN
    f15dat.nscreen = 0;
    % IHOT (cold start)
    f15dat.ihot = 0;
    % ICS 
    if all(abs(obj.p(:,1))) <= 180 && all(abs(obj.p(:,2))) <= 90
        f15dat.ics = 2 ; % uses lat, lon (probably)
    else
        f15dat.ics = 1 ; % not lat, lon
    end
    % IM (default is explicit with smagorinsky)
    f15dat.im = 511112 ;
    f15dat.iden = [];
    % NOLIBF (normal bottom friction)
    f15dat.nolibf = 1 ;
    % NOLIFA (fininte amplitude terms and wetting-drying on)
    f15dat.nolifa = 2 ;
    % NOLICA (the advection terms on)
    f15dat.nolica = 1 ; 
    % NOLICAT (the advection terms on)
    f15dat.nolicat = 1 ; 
    % NWP
    f15dat.nwp = 0 ;
    % NCOR 
    if f15dat.ics == 2
        f15dat.ncor = 1 ; %spatially varying coriolis from latitude
    else
        f15dat.ncor = 0 ; % read in Coriolis
    end
    % NTIP
    f15dat.ntip = 0 ;
    % NWS
    f15dat.nws = 0 ;
    f15dat.wtimnc = 0;
    % NRAMP
    f15dat.nramp = 0 ;
    % G
    f15dat.gravity = 9.81 ; 
    % TAU0 (the automatic spatially varying tau0)
    f15dat.tau0 = -1 ; 
    f15dat.tau0minmax = 0;
    % DTDP 
    f15dat.dtdp = 0 ; 
    % STATIM
    f15dat.statim = 0 ; 
    % REFTIM
    f15dat.reftim = 0 ; 
    % RNDY
    f15dat.rndy = 0;
    % DRAMP
    f15dat.dramp = 0;
    % A00, B00, C00
    f15dat.a00b00c00 = [0 1 0] ; 
    % H0
    f15dat.h0 = [0.1 0 0 0.01];
    % SLAM0, SFEA0
    f15dat.slam = mean(obj.p) ; 
    % CF
    f15dat.taucf = 0.0025 ;
    % ESLM (smag 0.05 coefficient)
    f15dat.elsm = -0.05;
    % CORI
    f15dat.cori = 0; 
    % NTIF
    f15dat.ntif = 0;
    % NBFR
    f15dat.nbfr = 0 ; 
    %  ANGINN
    f15dat.anginn = 90 ;
    % Normal flow external boundary condition (needs to be fixed)
    ibtype = [2 12 22 32 52] ;
    sm = 0 ;
    if ~isempty(obj.bd)
        for k = 1: length(ibtype)
            sm = sm + sum(~(obj.bd.ibtype - ibtype(k))) ;
        end
    end
    if ( sm > 0 ) 
        f15dat.nffr = 0 ;
    end
    % NOUTE, TOUTSE, TOUTFE, NSPOOLE
    f15dat.oute = zeros(1,4) ; 
    f15dat.elvstaloc = [];
    % NSTAE
    f15dat.nstae = 0 ;
    % NOUTV, TOUTV, TOUTFV, NSPOOLV
    f15dat.outv = zeros(1,4);
    f15dat.velstaloc = []; 
    % NSTAV
    f15dat.nstav = 0 ;
    % NSTAM, TOUSM, TOUFM, NSPOOLM 
    f15dat.outm = zeros(1,4) ; 
    f15dat.metstaloc = [];
    % NSTAM 
    f15dat.nstam = 0 ;
     
    % NOUTGE
    f15dat.outge = zeros(1,4) ;
    % NOUTGV
    f15dat.outgv = zeros(1,4) ;
    % NOUTGM 
    f15dat.outgm = zeros(1,4); 
    % NFREQ 
    f15dat.nfreq = 0 ;
    % THAS, THAF, NHAINC, FMV
    f15dat.outhar = zeros(1,4) ;
    % NHASE, NHASV, NHAGE, NHAGV
    f15dat.outhar_flag = zeros(1,4) ;
    % NHSTAR, NHSINC
    f15dat.nhstar = zeros(1,2) ;
    % ITITER, ISLDIA, CONVCR, ITMAX
    f15dat.ititer = [1 0 1e-10 25] ; 
    % Extra lines for NETCDF
    f15dat.nextraline = 10; 
    f15dat.extraline(1).msg = obj.title ;
    f15dat.extraline(2).msg = 'Notre Dame CHL';
    f15dat.extraline(3).msg = 'OceanMesh2D';
    f15dat.extraline(4).msg = 'History: None';
    f15dat.extraline(5).msg = 'https://github.com/CHLNDDEV/OceanMesh2D/';
    f15dat.extraline(6).msg = 'Comments: None' ;
    f15dat.extraline(7).msg = 'Host: Name';
    f15dat.extraline(8).msg = 'Metric, MSL';
    f15dat.extraline(9).msg = 'name@instit.edu';
    f15dat.extraline(10).msg = '';
    
    % Put into the msh class
    obj.f15 = f15dat;
end

%% Changing following parameters based on standard require inputs and 
%% populated properties of the msh class
% NSCREEN (once per simulation hour)
obj.f15.nscreen = round(3600/dt);
% DTDP (time step)
obj.f15.dtdp = dt ; 
% Runday
obj.f15.rndy = datenum(te)-datenum(ts);

obj.f15.extraline(10).msg = ts;

% Checking f11
if ~isempty(obj.f11)
    % using WJP's new IM for a 2DDI Baroclinic run
    obj.f15.im = 511113 ;
    obj.f15.iden = 1;
end 

% Checking f13
if ~isempty(obj.f13)
    obj.f15.nwp = obj.f13.nAttr ;
    if ( obj.f15.nwp > 0 ) 
        for l = 1: obj.f15.nwp
            obj.f15.AttrName(l).name = obj.f13.defval.Atr(l).AttrName ; 
            if strfind(obj.f15.AttrName(l).name,'primitive')
                obj.f15.tau0 = -3; % change tau0 option
            end
        end
    end
end

%% Test optional arguments
% default
const = [];
tidal_database = [];
sta_database = [];
if ~isempty(varargin)
    names = {'const','tidal_database','NWS','WTIMNC','tau0','tau0minmax',...
             'sta_database'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                const = varargin{ind*2}; 
            elseif ii == 2
                tidal_database = varargin{ind*2};
            elseif ii == 3
                obj.f15.nws = varargin{ind*2};
            elseif ii == 4
                obj.f15.wtimnc = varargin{ind*2};
            elseif ii == 5
                obj.f15.tau0 = varargin{ind*2};
            elseif ii == 6
                obj.f15.tau0minmax = varargin{ind*2};
            elseif ii == 7
                sta_database = varargin{ind*2};
            end
        end    
    end
end   

%% Changing following parameters based on optional argument inputs
% WTIMINC
if obj.f15.nws > 0 %~mod(obj.f15.nws - 1,10) ||  ~mod(obj.f15.nws - 11,10) || obj.f15.nws == 8 
    if obj.f15.wtimnc == 0
        error(['Must specify WTIMNC option with NWS = ' num2str(obj.f15.nws)])
    end
end

% Tau0FullDomainMin, Tau0FullDomainMax 
if ( abs(obj.f15.tau0 + 5.D0) < 1e-10 )
    if obj.f15.tau0minmax == 0
        error(['Must specify tau0minmax with tau0 option = ' num2str(obj.f15.tau0)])
    end
end

% Tidal potential and harmonic analysis stuff
if ~isempty(const)
    % Get the tidal factors etc based on ts and te 
    obj = tide_fac(obj,ts,te,const);
    
    % Harmonic analysis stuff (copy in tidal potential)
    obj.f15.nfreq = obj.f15.ntif; K = [];
    for k = 1: obj.f15.nfreq
        obj.f15.harfreq(k).name = obj.f15.tipotag(k).name;
        % get only the frequencies, nodal factor and equilibrium argument
        obj.f15.harfreq(k).val = obj.f15.tipotag(k).val([2,4,5]);
        if isnan(obj.f15.tipotag(k).val(1))
            % get rid of this constituent from potential
            obj.f15.ntif = obj.f15.ntif - 1;
            K(end+1) = k;
        end
    end
    obj.f15.tipotag(K) = [];
    
    if obj.f15.ntip == 0
        disp('Setting ntip = 1')
        obj.f15.ntip = 1;
    end
end

% Elevation Specified Boundary Conditions (tides)
if ~isempty(const) && ~isempty(tidal_database)
    % copy in the harmonic analysis stuff
    obj.f15.nbfr = obj.f15.nfreq;
    for k = 1: obj.f15.nbfr
        obj.f15.bountag(k).name = obj.f15.harfreq(k).name;
        % get only the frequencies, nodal factor and equilibrium argument
        obj.f15.bountag(k).val = obj.f15.harfreq(k).val;
    end

    % Do the interpolation to the boundaries
    obj = tidal_data_to_ob(obj,tidal_database,{obj.f15.bountag.name});
end

if ~isempty(const)
   % Here for the SA tidal potential the equilibrium argument is now
   % corrected to be referenced to the vernal equinox (subtract of about
   % 77 deg). But the equilibrim argument for the boundary and the analysis
   % is kept as referenced to the Calender year.
   ii = find(startsWith({obj.f15.tipotag.name},'SA'));
   if ~isempty(ii)
       % need to refer to vernal equinox instead of calender year, 
       % at least for the tidal potential (subtract about 76.9 deg);
       obj.f15.tipotag(ii).val(5) = obj.f15.tipotag(ii).val(5) - 78/365*360;
       if obj.f15.tipotag(ii).val(5) < 0
          obj.f15.tipotag(ii).val(5) = obj.f15.tipotag(ii).val(5) + 360; 
       end
   end
end
% 
if ~isempty(sta_database)
    % parse the station reader
    obj = tidal_stations_parser(obj,string(sta_database(1:2:end-1)),...
                                sta_database(2:2:end));
end
% 
% % NOUTGC
% if ( f15dat.im == 10 )
%     f15dat.outgc = readlinevec( fid ) ; 
% end
% 
% % NOUTGW
% if ( f15dat.nws ~= 0 ) 
%     f15dat.outgw = readlinevec( fid ) ;
% end

end