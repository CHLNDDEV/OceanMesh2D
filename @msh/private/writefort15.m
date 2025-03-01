function fid = writefort15( f15dat, f15out, boudat )
%
if ( strcmp(strtrim(f15out), 'fort.15') )
    disp('Error: an output file name must not be fort.15') ;
    return ;
end

fid = fopen( f15out, 'w' ) ;

fprintf(fid, '%s\n', f15dat.rundes ) ; % RUNDES
fprintf(fid, '%s\n', f15dat.runid  ) ;  % RUNID

% NFOVER
L = length(f15dat.nfover);
dfmt = repmat('%g ',1, L);
fprintf( fid, [dfmt '        \t ! NFOVER \n'], f15dat.nfover ) ;

% NABOUT
fprintf(fid, '%d        \t ! NABOUT \n', f15dat.nabout ) ;

% NSCREEN
fprintf(fid, '%d        \t ! NSCREEN \n', f15dat.nscreen ) ;

% IHOT
fprintf(fid, '%d        \t ! IHOT \n', f15dat.ihot ) ;

% ICS
fprintf(fid, '%d        \t ! ICS \n', f15dat.ics ) ;

% IM
fprintf(fid, '%d        \t ! IM \n', f15dat.im ) ;

% IDEN
if f15dat.im == 20 || f15dat.im == 30
    fprintf( fid, '%d           \t ! IDEN \n', f15dat.iden ) ;
end

% NOLIBF
fprintf(fid, '%d        \t ! NOLIBF \n', f15dat.nolibf ) ;

% NOLIFA
fprintf(fid, '%d        \t ! NOLIFA \n', f15dat.nolifa ) ;

% NOLICA
fprintf(fid, '%d        \t ! NOLICA \n', f15dat.nolica ) ;

% NOLICAT
fprintf(fid, '%d        \t ! NOLICAT \n', f15dat.nolicat ) ;

% NWP
fprintf(fid, '%d        \t ! NWP \n', f15dat.nwp ) ;
if ( f15dat.nwp > 0 )
    for l = 1: f15dat.nwp
        fprintf(fid, '%s\n', f15dat.AttrName(l).name ) ;
    end
end

% NCOR
fprintf(fid, '%d   \t \t ! NCOR \n', f15dat.ncor ) ;

% NTIP
fprintf(fid, '%d   \t \t ! NTIP \n', f15dat.ntip ) ;

% NWS
fprintf(fid, '%d   \t \t ! NWS \n', f15dat.nws ) ;

% NRAMP
fprintf(fid, '%d   \t \t ! NRAMP \n', f15dat.nramp ) ;

% G
fprintf(fid, '%f   \t ! G  \n', f15dat.gravity ) ;

% TAU0
fprintf(fid, '%g   \t \t ! TAU0 \n', f15dat.tau0 ) ;

% Tau0FullDomainMin, Tau0FullDomainMax
if ( abs(f15dat.tau0 + 5.D0) < 1e-10 )
    fprintf(fid, '%f %f   \t ! Tau0FullDomainMin, Tau0FullDomainMax \n', f15dat.tau0minmax ) ;
end

% DTDP
fprintf( fid, '%g  \t \t ! DTDP \n', f15dat.dtdp ) ;

% STATIM
fprintf( fid, '%g  \t \t ! STATIM \n',  f15dat.statim ) ;

% REFTIM
fprintf( fid, '%g  \t \t ! REFTIM \n', f15dat.reftim ) ;

% WTIMINC
if f15dat.nws > 0
    nws = f15dat.nws;
    wtimnc = f15dat.wtimnc;
    rstimnc_fmt = ''; rstimnc_msg = '';
    if nws > 300
        wtimnc = [f15dat.wtimnc f15dat.rstimnc];
        nws = f15dat.nws - 300;
        rstimnc_fmt = ' %d';
        rstimnc_msg = ' RSTIMINC';
    end
    if nws == 8
        wtimnc_fmt = '%d %d %d %d %d %g';
        wtimnc_msg = 'YYYY MM DD HH24 StormNumber BLAdj';
    elseif nws >= 19
        wtimnc_fmt = '%d %d %d %d %d %g %d';
        wtimnc_msg = 'YYYY MM DD HH24 StormNumber BLAdj geofactor';
    else
        wtimnc_fmt = '%d  ';
        wtimnc_msg = 'WTIMINC';
    end
    wtimnc_fmt = [wtimnc_fmt rstimnc_fmt];
    wtimnc_msg = [wtimnc_msg rstimnc_msg];
    fprintf( fid, wtimnc_fmt, wtimnc ) ;
    fprintf( fid, ['  \t ! ' wtimnc_msg ' \n'] ) ;
end

% RNDY
fprintf( fid, '%g   \t \t ! RNDY \n', f15dat.rndy ) ;

% DRAMP
L = length(f15dat.dramp);
dfmt = repmat('%g ',1, L);
fprintf( fid, [dfmt ' \t \t ! DRAMP \n'],  f15dat.dramp ) ;

% A00, B00, C00
fprintf( fid, '%g %g %g  \t ! A00, B00, C00 \n', f15dat.a00b00c00 ) ;

% H0
len = length(f15dat.h0) ;
for k = 1: len
    fprintf( fid, '%g ', f15dat.h0(k) ) ;
end
fprintf( fid, '   \t ! H0, 2*dummy, VELMIN \n' ) ;

% SLAM0, SFEA0
fprintf( fid, '%f %f  \t \t ! SLAM0, SFEA0 \n', f15dat.slam ) ;

% CF
if ( f15dat.nolibf <= 2 )
    fprintf( fid, '%f ', f15dat.taucf ) ;
    fprintf( fid, '    \t ! CF \n' ) ;
end

% ESLM, ESLC
if ( f15dat.im <= 2 || f15dat.im == 10 || f15dat.im >= 111111 )
    fprintf( fid, '%f ', f15dat.elsm ) ;
    fprintf( fid, '    \t ! ELSM \n' ) ;
end

% CORI
fprintf( fid, '%f   \t ! CORI \n', f15dat.cori ) ;

% NTIF
fprintf( fid, '%d   \t \t ! NTIF \n', f15dat.ntif ) ;

% Tidal potential
for k = 1: f15dat.ntif
    fprintf( fid, '%s \n',  f15dat.tipotag(k).name ) ;
    fprintf( fid, '%f %16.9e %f %f %f', f15dat.tipotag(k).val ) ;
    fprintf( fid, '   \t !  TPK, AMIGT, ETRF, FFT, FACET \n' ) ;
end

% NBFR
fprintf( fid, '%d   \t \t ! NBFR \n', f15dat.nbfr ) ;
for k = 1: f15dat.nbfr
    fprintf( fid, '%s  \n', f15dat.bountag(k).name ) ;
    fprintf( fid, '%16.9e %f %f \n', f15dat.bountag(k).val ) ;
end

% Open boundary harmonic forcing
for k = 1: f15dat.nbfr
    fprintf(fid, '%s \n', f15dat.opealpha(k).name  ) ;
    fprintf(fid, '%16.9e %16.10g \n', f15dat.opealpha(k).val' ) ;
end

%  ANGINN
fprintf( fid, '%g   \t \t ! ANGINN \n', f15dat.anginn ) ;

% Land boundary
ibtype = [2 12 22 32 52] ;
sm = 0 ;

if ~isempty(boudat)
    for k = 1: length(ibtype)
        sm = sm + sum(~(boudat.ibtype - ibtype(k))) ;
    end
end
if ( sm > 0 )
    fprintf( fid, '%d  \t \t ! NFFR \n', f15dat.nffr ) ;
    
    nm = 0 ;
    for ib = 1: boudat.nbou
        ibty = boudat.ibtype(ib);
        
        switch ibty
            case {2,12,22,32,52}
                nm = nm + boudat.nvell(ib) ;
            otherwise
        end
    end
    
    for k = 1: f15dat.nffr
        fprintf( fid, '%s \n', f15dat.fbountag(k).name ) ;
        fprintf( fid, '%15.8e ', f15dat.fbounspec(k).val ) ;
        fprintf( fid, '\n') ;
    end
    
    
    for k = 1: f15dat.nffr
        fprintf( fid, '%s  \n',  f15dat.boualpha(k).name  ) ;
        
        % val = fscanf(fid, '%f ' ) ; % Must be revisit
        icnt = 0 ;
        for ib = 1: boudat.nbou
            ibty = boudat.ibtype(ib);
            
            switch ibty
                case {2,12,22,52}
                    for ir = 1: boudat.nvell(ib)
                        icnt = icnt + 1 ;
                        
                        fprintf(fid, '%16.9e  ', f15dat.qnam(k).val(icnt,1:2) ) ;
                        fprintf(fid, ' \n' ) ;
                    end
                case 32
                    for ir = 1: boudat.nvell(ib)
                        icnt = icnt + 1 ;
                        
                        fprintf(fid, '%15.8e ', f15dat.qnam(k).val(icnt,1:5) ) ;
                        fprintf(fid, ' \n' ) ;
                    end
                otherwise
            end
        end
    end
end

% NOUTE, TOUTSE, TOUTFE, NSPOOLE
fprintf( fid, '%d %g %g %d', f15dat.oute ) ;
fprintf( fid, '   \t ! NOUTE, TOUTSE, TOUTFE, NSPOOLE \n' ) ;

% NSTAE
fprintf( fid, '%d  \t \t ! NSTAE \n', f15dat.nstae ) ;

% STAE location
if ( f15dat.nstae > 0 )
    for k = 1: f15dat.nstae
        fprintf(fid, '%f %f %s \n', f15dat.elvstaloc(k,1:2), f15dat.elvstaname{k} ) ;
    end
end

% NOUTV, TOUTV, TOUTFV, NSPOOLV
fprintf( fid, '%d %g %g %d',  f15dat.outv ) ;
fprintf( fid, '  \t ! NOUTV, TOUTV, TOUTFV, NSPOOLV \n') ;

% NSTAV
fprintf( fid, '%d  \t \t ! NSTAV \n', f15dat.nstav ) ;

% STAV location
if ( f15dat.nstav > 0 )
    for k = 1: f15dat.nstav
        fprintf( fid, '%f %f %s \n', f15dat.velstaloc(k,1:2), f15dat.velstaname{k} ) ;
    end
end

% NOUTM, TOUTM, TOUTFM, NSPOOLM
if ( f15dat.nws ~= 0  )
    fprintf(fid, '%d %g %g %d', f15dat.outm ) ;
    fprintf(fid, '  \t ! NOUTM, TOUTM, TOUTFM, NSPOOLM \n') ;
    
    % NSTAM
    fprintf(fid, '%d   \t \t ! NSTAM \n', f15dat.nstam ) ;
    
    % STAM location
    if ( f15dat.nstam > 0 )
        for k = 1: f15dat.nstam
            fprintf( fid, '%f %f %s \n', f15dat.metstaloc(k,1:2), f15dat.metstaname{k} ) ;
        end
    end
end

% NOUTGE
fprintf( fid, '%d %g %g %d', f15dat.outge ) ;
fprintf( fid, '  \t ! NOUTGE, ... \n') ;

% NOUTGV
fprintf( fid, '%d %g %g %d', f15dat.outgv ) ;
fprintf( fid, '  \t ! NOUTGV, ... \n') ;

% NOUTGC
if ( f15dat.im == 10 )
    fprintf( fid, '%d %g %g %d', f15dat.outgc ) ;
    fprintf( fid, '  \t ! NOUTGC, ... \n' ) ;
end

% NOUTGM
if ( f15dat.nws ~= 0  )
    fprintf( fid, '%d %g %g %d', f15dat.outgm ) ;
    fprintf( fid, '  \t ! NOUTGM, ... \n' ) ;
end

% NFREQ
fprintf( fid, '%d  \t \t ! NFREQ \n', f15dat.nfreq ) ;
for k = 1: f15dat.nfreq
    fprintf( fid, '%s \n', f15dat.harfreq(k).name ) ;
    
    fprintf( fid, '%16.9e %f %f \n', f15dat.harfreq(k).val ) ;
end

% THAS, THAF, NHAINC, FMV
fprintf( fid, '%g %g %d %g  \t ! THAS, THAF, NHAINC, FMV \n', f15dat.outhar ) ;

% NHASE, NHASV, NHAGE, NHAGV
fprintf( fid, '%d ', f15dat.outhar_flag ) ;
fprintf( fid, '  \t ! NHASE, NHASV, NHAGE, NHAGV \n' ) ;

% NHSTAR, NHSINC
fprintf( fid, '%d %d  \t ! NHSTAR, NHSINC \n', f15dat.nhstar) ;

% ITITER, ISLDIA, CONVCR, ITMAX
fprintf( fid, '%d %d  %16.9e %d  \t ! ITITER, ISLDIA, CONVCR, ITMAX \n', f15dat.ititer ) ;

% Extra lines including NETCDF & namelist
for k = 1: f15dat.nextraline
    fprintf( fid, '%s\n', f15dat.extraline(k).msg ) ;
end

% If
if find(strcmp(fieldnames(f15dat),'controllist'),1)

for k = 1: length(f15dat.controllist)
    fprintf( fid, '! -- Begin %s Control Namelist -- \n', f15dat.controllist(k).type ) ;
    fprintf( fid, '&%sControl\n', f15dat.controllist(k).type ) ;
    for m = 1:length(f15dat.controllist(k).var)
        val = f15dat.controllist(k).var(m).val;
        if ~ischar(val); 
           if islogical(val)
              if val; val = 'T'; else; val = 'F'; end
           else
              val = num2str(val);
           end
           fprintf( fid, '%s = %s,\n',f15dat.controllist(k).var(m).name,val);
        else
           fprintf( fid, '%s = "%s",\n',f15dat.controllist(k).var(m).name,val);
        end
    end
    fprintf( fid, '/ ! -- End %s Control Namelist -- \n', f15dat.controllist(k).type ) ;
end

fclose(fid) ;

end
