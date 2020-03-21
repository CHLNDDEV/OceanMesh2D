function Error = Calc_TideGaugeErrors(fort53,stadat,varargin)
% Calculates the complex root mean square error of the major 8
% constitutents given a station database and writes a text file called
% Harm.txt. This file is then used with HarmonicsCompare to validate the
% tides by producing scatter plots and xml files for displaying on a web
% server with google maps. 

% authors: william pringle and keith roberts. 
% CHL, UND, March 2018.

%%
if nargin < 2, error('Not enough input arguments!'); end
const_cell = {'M2','S2','N2','K2','K1','O1','P1','Q1','all'};
type_cell = {'Pelagic','Shelf','coast'};
write_out = true;
if nargin > 2
    const_cell = varargin{1};
end
if nargin > 3
   type_cell = varargin{2};
end
if nargin > 4
    write_out = varargin{3};
end

%% 
x = ncread(fort53,'x');
y = ncread(fort53,'y');
depth = ncread(fort53,'depth');

mdl2 = KDTreeSearcher([x,y]);
x = x'; y = y';

c_list = {'M2','S2','N2','K2','K1','O1','P1','Q1'};
c_list2 = ncread(fort53,'const');

tt = 0;
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
        RMS = zeros(length(I),1); V_o = 0*RMS; V_m = 0*RMS;
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
            V_o = V_o + a_0.^2;
            V_m = V_m + a_m.^2;
        end
        K = find(RMS == 0);
        latt = T_o.Lat; lont = T_o.Lon;
        a_m(K) = []; a_0(K) = [];
        g_m(K) = []; g_0(K) = [];
        latt(K) = []; lont(K) = [];
        RMS(K) = []; V_m(K) = []; V_o(K) = [];
        %RSS = sqrt(0.5*sum(RMS)/length(RMS))*100
        RMS = sqrt(0.5*RMS);
        V_o = sqrt(0.5*V_o);
        V_m = sqrt(0.5*V_m);
        
        Error(tt,cc).Const = const{1};
        Error(tt,cc).Type = type{1};
        Error(tt,cc).RMSE = RMS;
        Error(tt,cc).V_o = V_o;
        Error(tt,cc).V_m = V_m;
        Error(tt,cc).Lat = latt;
        Error(tt,cc).Lon = lont;  
        if ~strcmp(const{1},'all')
            Error(tt,cc).A_obs = a_0;
            Error(tt,cc).A_mod = a_m;    
            Error(tt,cc).G_obs = g_0;
            Error(tt,cc).G_mod = g_m;       
        end
    end
end

if ~write_out; return; end

% Create Harm.txt file for Harmonic_compare to work properly
fid=fopen('Harm.txt','w');  %<--this will be what this script produces!!
CLASS = length(type);   % class == 1 deep/pelagic, class == 2 shelf, class == 3 nearshore
CONST = [1:length(const_cell)]; % index into constitutent const_cell array
NoSta = length(Error(CLASS,CONST(1)).RMSE);
NoConst= length(CONST);
%fullmarkers={'s','o','h','^','v','x','+','d'};

fprintf(fid,'%d %d\n',NoSta,NoConst) ; %<--change number of stations for type of stations
fprintf(fid,'%s %s %s %s %s %s %s %s\n','s','o','h','^','v','x','+','d');  %<-- depending on length of const_cell, add more marker types.

for sta = 1 : NoSta
    %fprintf(fid,'%s %f %f\n', char(Sta_Names{CLASS,CONST(1)}(sta)),Lon_cell{CLASS,CONST(1)}(sta),Lat_cell{CLASS,CONST(1)}(sta));
    fprintf(fid,'%s %f %f\n', num2str(sta),Error(CLASS,CONST(1)).Lon(sta),Error(CLASS,CONST(1)).Lat(sta));
    
    for c = 1 : NoConst
        constituent = const_cell{CONST(c)};
        if strcmp(const,'all'); continue; end
        %   const-name   obs-amp     mod-amp     obs-phs     mod-phs
        fprintf(fid,'%s %f %f %f %f\n',constituent,Error(CLASS,CONST(c)).A_obs(sta),Error(CLASS,CONST(c)).A_mod(sta),...
            Error(CLASS,CONST(c)).g_obs(sta),Error(CLASS,CONST(c)).g_mod(sta));
    end
end
fclose(fid);
%EOF
end
