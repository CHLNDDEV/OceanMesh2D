function f15dat = readfort15( f15name, opedat, boudat )
%
%
% Readfort.15
%    
%

f15file = f15name ; 

fid = fopen( strtrim(f15file) ) ;

f15dat.rundes = fgetl(fid) ; % RUNDES
f15dat.runid = fgetl(fid) ;  % RUNID

% NFOVER
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.nfover = str2num(token) ;

% NABOUT
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.nabout = str2num(token) ;

% NSCREEN
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.nscreen = str2num(token) ;

% IHOT
line = fgetl(fid) ;
[token,res] = strtok(line) ;
f15dat.ihot = str2num(token) ;

% ICS
line = fgetl(fid) ;
[token,res] = strtok(line) ;
f15dat.ics = str2num(token) ;

% IM
f15dat.iden = [] ;
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.im = str2num(token) ;

% IDEN
if ( f15dat.im == 20 || f15dat.im == 30 || f15dat.im == 511113 ) 
    line = fgetl(fid) ;
    [token,res] = strtok(line) ; 
    
    f15dat.iden = str2num(token) ;  
end

% NOLIBF
token = readlinetoken( fid ) ;
f15dat.nolibf = str2num(token) ;

% NOLIFA
token = readlinetoken( fid ) ;
f15dat.nolifa = str2num(token) ;

% NOLICA
token = readlinetoken( fid ) ;
f15dat.nolica = str2num(token) ; 

% NOLICAT
token = readlinetoken( fid ) ; 
f15dat.nolicat = str2num( token ) ; 

% NWP
f15dat.nwp = 0 ;
token = readlinetoken( fid ) ;
f15dat.nwp = str2num( token ) ;
if ( f15dat.nwp > 0 ) 
    for l = 1: f15dat.nwp
        f15dat.AttrName(l).name = fgetl(fid) ; 
    end
end

% NCOR
token = readlinetoken( fid ) ;
f15dat.ncor = str2num( token ) ; 

% NTIP
token = readlinetoken( fid ) ; 
f15dat.ntip = str2num( token ) ;

% NWS
token = readlinetoken( fid ) ; 
f15dat.nws = str2num( token ) ;

% NRAMP
token = readlinetoken( fid ) ; 
f15dat.nramp = str2num( token ) ; 

% G
token = readlinetoken( fid ) ; 
f15dat.gravity = str2num( token ) ; 

% TAU0
token = readlinetoken( fid ) ;
f15dat.tau0 = str2num( token ) ; 

% Tau0FullDomainMin, Tau0FullDomainMax 
if ( abs(f15dat.tau0 + 5.D0) < 1e-10 )
    f15dat.tau0minmax =  readlinevec( fid ) ;
end

% DTDP
token = readlinetoken( fid ) ; 
f15dat.dtdp = str2num( token ) ; 

% STATIM
token = readlinetoken( fid ) ; 
f15dat.statim = str2num( token ) ; 

% REFTIM
token = readlinetoken( fid ) ; 
f15dat.reftim = str2num( token ) ; 

% WTIMINC
if f15dat.nws > 0
   f15dat.wtimnc = readlinevec( fid ) ;
end
  
% RNDY
token = readlinetoken( fid ) ; 
f15dat.rndy = str2num( token ) ;

% DRAMP
f15dat.dramp = readlinevec( fid ) ;

% A00, B00, C00
f15dat.a00b00c00 = readlinevec( fid ) ; 

% H0
f15dat.h0 = readlinevec( fid ) ;

% SLAM0, SFEA0
f15dat.slam = readlinevec( fid ) ; 

% TAU, CF, HBREAK, FTHETA, FGAMMA
if ( f15dat.nolibf <= 2 ) 
    f15dat.taucf = readlinevec( fid ) ;
end

% ESLM, ESLC
f15dat.elsm = [] ;
if ( f15dat.im <= 2 | f15dat.im == 10 | f15dat.im >= 111111 )
    f15dat.elsm = readlinevec( fid ) ; 
end

% CORI
f15dat.cori = readlinescalar( fid ) ; 

% NTIF
f15dat.ntif = readlinescalar( fid ) ;

% Tidal potential
for k = 1: f15dat.ntif
    f15dat.tipotag(k).name = fgetl( fid ) ;
    f15dat.tipspec(k).val =  readlinevec( fid ) ; 
end
 
% NBFR
f15dat.nbfr = readlinescalar( fid ) ; 
for k = 1: f15dat.nbfr
    f15dat.bountag(k).name = fgetl( fid ) ; 
    f15dat.bounspec(k).val = readlinevec( fid ) ; 
end

% Open boudnary harmonic forcing  
for k = 1: f15dat.nbfr
    f15dat.opealpha(k).name = fgetl( fid ) ; 
    
    nvd = opedat.neta ; 
    
    val = fscanf(fid, '%f %f \n', [2 nvd] ) ; 
    
    %len = size(val) 
    % val = reshape( val, 2, len/2 ) ; 
    f15dat.opeemoefa(k).val = val' ; 
end


%  ANGINN
f15dat.anginn = readlinescalar( fid ) ;

% Land boundary 
ibtype = [2 12 22 32 52] ;
sm = 0 ;

for k = 1: length(ibtype)
    sm = sm + sum(~(boudat.ibtype - ibtype(k))) ;
end
if ( sm > 0 ) 
    f15dat.nffr = readlinescalar( fid ) ;
    
    nm = 0 ; 
    for ib = 1: boudat.nbou
        ibty = boudat.ibtype(ib)
        
         switch ibty
             case {2,12,22,52}
                 nm = nm + boudat.nvell(ib) ;
             case 32
                 nm = nm + boudat.nvell(ib) ;
             otherwise
         end
    end
    
    for k = 1: f15dat.nffr
        f15dat.fbountag(k).name = fgetl( fid ) ; 
        f15dat.fbounspec(k).val = readlinevec( fid ) ; 
    end
    
    for k = 1: f15dat.nffr
        f15dat.qnam(k).val = zeros(nm,5) ; 
        f15dat.boualpha(k).name = fgetl( fid ) ; 
        
        % val = fscanf(fid, '%f ' ) ; % Must be revisit
        icnt = 0 ;
        for ib = 1: boudat.nbou
            ibty = boudat.ibtype(ib) ;
            
            switch ibty
                case {2,12,22,52}
                   for ir = 1: boudat.nvell(ib)
                      val = readlinevec( fid ) ;
                  
                      icnt = icnt + 1 ;
                      f15dat.qnam(k).val(icnt,1:2) = val ; 
                   end
                case 32
                   for ir = 1: boudat.nvell(ib)
                      val = readlinesvec( fid ) ;
                  
                      icnt = icnt + 1 ;
                      f15dat.qnam(k).val(icnt,1:5) = val ; 
                   end 
                otherwise
            end
        end
    end
end

% NOUTE, TOUTSE, TOUTFE, NSPOOLE
f15dat.oute = readlinevec( fid ) ; 

% NSTAE
f15dat.nstae = readlinescalar( fid ) ;
 
% STAE location
if ( f15dat.nstae > 0 )
    
    f15dat.elvstaloc = zeros(f15dat.nstae,2) ;
    for k = 1: f15dat.nstae
        val = readlinevec( fid ) ;
        
        f15dat.elvstaloc(k,1:2) = val(1:2)' ;
    end
end

% NOUTV, TOUTV, TOUTFV, NSPOOLV
f15dat.outv = readlinevec( fid ) ;

% NSTAV
f15dat.nstav = readlinescalar( fid ) ;
f15dat.velstaloc = [] ;

% STAV location
if ( f15dat.nstav > 0 ) 
    f15dat.velstaloc = zeros(f15dat.nstav,2) ;
    for k = 1: f15dat.nstav
        val = readlinevec( fid ) ; 
        
        f15dat.velstaloc(k,1:2) = val ; 
    end
end

% NOUTM, TOUTM, TOUTFM, NSPOOLM
if ( f15dat.nws ~= 0 )
    f15dat.outm = readlinevec( fid ) ;
    
    % NSTAM
    f15dat.nstam = readlinescalar( fid ) ; 
    f15dat.metstaloc = [] ;
    
    % STAM location
    if ( f15dat.nstam > 0 )
        f15dat.metstaloc = zeros(f15dat.nstam,2) ;
        for k = 1: f15dat.nstam
            val = readlinevec( fid ) ;
            
            f15dat.metstaloc(k,1:2) = val ;
        end
    end
end

% NOUTGE
f15dat.outge = readlinevec( fid ) ;

% NOUTGV
f15dat.outgv = readlinevec( fid ) ;

% NOUTGC
if ( f15dat.im == 10 )
    f15dat.outgc = readlinevec( fid ) ; 
end

% NOUTGW
if ( f15dat.nws ~= 0 ) 
    f15dat.outgw = readlinevec( fid ) ;
end


% NFREQ 
f15dat.nfreq = readlinescalar( fid ) ; 
for k = 1: f15dat.nfreq
    f15dat.harfreq(k).name = fgetl( fid ) ; 
    
    f15dat.harfreq(k).val =  readlinevec( fid ) ;
end

% THAS, THAF, NHAINC, FMV
f15dat.outhar = readlinevec( fid ) ;

% NHASE, NHASV, NHAGE, NHAGV
f15dat.outhar_flag = readlinevec( fid ) ;

% NHSTAR, NHSINC
f15dat.nhstar = readlinevec( fid ) ;

% ITITER, ISLDIA, CONVCR, ITMAX
f15dat.ititer = readlinevec( fid ) ; 

%
% Extra liens including NETCDF & namelist 
%
line = fgetl( fid ) ;
icount = 0 ;
while ischar(line)
    icount = icount + 1 ;
    f15dat.extraline(icount).msg = line ; 
    
    line = fgetl( fid ) ;
end
f15dat.nextraline = icount ; 


fclose(fid) ;

    function vec = readlinevec( fid )
        msg = fgetl(fid) ; 
        
        idx = find(msg == ',') ; 
        if ( ~isempty(msg) ) 
            msg(idx) = ' ' ;
        end
        vec = sscanf(msg,'%f') ; 
    end

    function val = readlinescalar( fid )
        
       msg = fgetl(fid) ; 
       
       [ftoken,rem] = strtok(msg) ; 
       
       val = str2num(ftoken) ;
    end

    function token = readlinetoken( fid )
        
       msg = fgetl(fid) ; 
       [token,rem] = strtok(msg) ; 
        
    end
 
end