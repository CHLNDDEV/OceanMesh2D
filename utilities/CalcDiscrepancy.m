function [fid]=ValidateTides(fort53,stadat,varargin)
% Calculates the complex root mean square error of the major 8
% constitutents given a station database and writes a text file called
% Harm.txt. This file is then used with HarmonicsCompare to validate the
% tides by producing scatter plots and xml files for displaying on a web
% server with google maps. 

% authors: william pringle and keith roberts. 
% CHL, UND, March 2018.

%%
if nargin < 2, error('Not enough input arguments!'); end;
if nargin ==3
   % user defines constitutents 
else
    const_cell = {'M2','S2','N2','K2','K1','O1','P1','Q1','all'};
end
if nargin == 4
  % user defines type_cell  
else
    type_cell = {'Pelagic','Shelf','coast'};
end

%%
x = ncread(fort53,'x');
y = ncread(fort53,'y');
depth = ncread(fort53,'depth');

mdl2 = KDTreeSearcher([x,y]);
x = x'; y = y';

c_list2 = ncread(fort53,'const');

tt = 0;
fid=fopen('Harm.txt','w');  %<--this will be what this script produces!!
for type = type_cell
    tt = tt + 1;
    T_o = readtable(stadat);
    if strcmp(type{1},'coast')
        I1 = strfind(T_o.Source,'Truth');
        I2 = strfind(T_o.Name,'_TP');
        I3 = strfind(T_o.Name,'_Shelf');
        I = find(cellfun(@isempty,I1) & cellfun(@isempty,I2) & cellfun(@isempty,I3));
    else
        I1 = strfind(T_o.Source,type{1});
        I2 = strfind(T_o.Name,type{1});
        I = find(~cellfun(@isempty,I1) | ~cellfun(@isempty,I2));
    end
    
    T_o = T_o(I,:);
    
    idx1 = knnsearch(mdl2,[T_o.Lon T_o.Lat],'k',12);
    cc = 0;
    for const = const_cell
        cc = cc + 1;
        pos = 0;
        RMS = zeros(length(I),1); V = zeros(length(I),1);
        for c = c_list
            pos = pos + 1;
            if ~strcmp(const{1},'all') && ~strcmp(c{1},const{1})
                continue
            end
            a_0 = table2array(T_o(:,4+2*pos-1));
            g_0 = table2array(T_o(:,4+2*pos));
            
            pos2 = 0;
            for ccc = c_list2
                pos2 = pos2 + 1;
                if strcmp(strtrim(ccc'),c{1})
                    break
                end
            end
            
            % Using fort.53
            amp1 = ncread(fort53,'amp',[pos2 1],[1 length(x)]);
            phs1 = ncread(fort53,'phs',[pos2 1],[1 length(x)]);
            
            % These are nodes to make the interpolation
            [a_m, g_m] = sta_interp(x(idx1),y(idx1),depth(idx1),amp1(idx1),...
                phs1(idx1),T_o.Lon,T_o.Lat);
            
            % Do not add to sum
            a_m(isnan(a_0)) = 0; g_m(isnan(a_0)) = 0;
            g_0(isnan(a_0)) = 0; a_0(isnan(a_0)) = 0;
            g_0(isnan(a_m)) = 0; a_0(isnan(a_m)) = 0;
            g_m(isnan(a_m)) = 0; a_m(isnan(a_m)) = 0;
            
            %RMS = RMS + (a_0.^2 + a_m.^2 - 2*a_0.*a_m.*cosd(g_0 - g_m));
            RMS = RMS + (a_0.^2 + a_m.^2 - 2*a_0.*a_m.*cosd(g_0 - g_m));
            V   = V   + a_0.^2;
        end
        K = find(RMS == 0);
        latt = T_o.Lat; lont = T_o.Lon;
        latt(K) = []; lont(K) = [];
        RMS(K) = []; V (K) = [];
        %RSS = sqrt(0.5*sum(RMS)/length(RMS))*100
        RMS = sqrt(0.5*RMS);
        V = sqrt(0.5*V);
        
        D(tt,cc) = 100*mean(RMS);
        SD(tt,cc) = 100*std(RMS);
        %disp(D*100)
        
        RD(tt,cc) = 100*mean(RMS./V);
        SRD(tt,cc) = 100*std(RMS./V);
        
        D_cell{tt,cc} = RMS;
        a_m(K) = []; a_0(K) = [];
        OBSERVED{tt,cc} = a_0;
        MODELED{tt,cc}  = a_m;
        V_cell{tt,cc} = V;
        Lat_cell{tt,cc} = latt;
        Lon_cell{tt,cc} = lont;
        
    end
end

% Create Harm.txt file for Harmonic_compare to work properly
CLASS = 3;   % class == 1 deep/pelagic, class == 2 shelf, class == 3 nearshore
CONST = [1:8]; % index into constitutent const_cell array
NoSta = length(D_cell{CLASS,CONST(1)});
NoConst= length(CONST);
%fullmarkers={'s','o','h','^','v','x','+','d'};

fprintf(fid,'%d %d\n',NoSta,NoConst) ; %<--change number of stations for type of stations
fprintf(fid,'%s %s %s %s %s %s %s %s\n','s','o','h','^','v','x','+','d');  %<-- depending on length of const_cell, add more marker types.

for sta = 1 : NoSta
    %fprintf(fid,'%s %f %f\n', char(Sta_Names{CLASS,CONST(1)}(sta)),Lon_cell{CLASS,CONST(1)}(sta),Lat_cell{CLASS,CONST(1)}(sta));
    fprintf(fid,'%s %f %f\n', num2str(sta),Lon_cell{CLASS,CONST(1)}(sta),Lat_cell{CLASS,CONST(1)}(sta));
    
    for c = 1 : NoConst
        constituent = const_cell{CONST(c)};
        
        %   const-name   obs-amp     mod-amp     obs-phs     mod-phs
        fprintf(fid,'%s %f %f %f %f\n',constituent,OBSERVED_AMP{CLASS,CONST(c)}(sta),MODELED_AMP{CLASS,CONST(c)}(sta),...
            OBSERVED_PHS{CLASS,CONST(c)}(sta),MODELED_PHS{CLASS,CONST(c)}(sta));
    end
end

end
