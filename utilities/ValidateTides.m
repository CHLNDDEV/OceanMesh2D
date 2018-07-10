function [error_matrix,error_values,sta_locs,OBSERVED_AMP,MODELED_AMP] = ...
                                      ValidateTides(fort53,stadat,varargin)
% [error_matrix,error_values,sta_locs,OBSERVED_AMP,MODELED_AMP] = 
%                                     ValidateTides(fort53,stadat,varargin)
% Calculates the complex root mean square error of the major 8
% constitutents given a station database and writes a text file called
% Harm.txt. This file is then used with HarmonicsCompare to validate the
% tides by producing scatter plots and xml files for displaying on a web
% server with google maps.
% INPUTS:
% fort53: name of the fort.53.nc file.
% stadat: name of the station database to validate against.
% const_cell (optional): cell array of constitutent names to validate.
%                         by default: {'M2','S2','N2','K2','K1','O1','P1','Q1','all'};
% NOTE that last entry can be "all" to denote the mean omplex RMSE for all stations over all constitutents. 
% type_cell (optional): cell array describing the station types (Pelagic,
%                       shelf or coast).
%                       by default: {'Pelagic','Shelf','coast'};
%
% OUTPUTS
% error matrix: cell array of the mean complex rmse in the order of type of station by constitutent. 
% error_values: cell array of the complex rmse at each station in the order of type of station by constitutent. 
% REQUIRES: stainterp.m
%
% authors: william pringle and keith roberts.
% CHL, UND, March 2018.

%%
if nargin < 2, error('Not enough input arguments!'); end
if nargin >2 
    % user defines constitutents
    const_cell = varargin{1}; 
else
    const_cell = {'M2','S2','N2','K2','K1','O1','P1','Q1'};
end
if nargin > 3
    % user defines type_cell
    type_cell = varargin{2}; 
else
    type_cell = {'Pelagic','Shelf','coast'};
end
c_list = {'M2','S2','N2','K2','K1','O1','P1','Q1'}; 
%%
x = ncread(fort53,'x');
y = ncread(fort53,'y');
depth = ncread(fort53,'depth');

mdl2 = KDTreeSearcher([x,y]);
x = x'; y = y';

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
        %a_m(K) = []; a_0(K) = [];
        
        OBSERVED_AMP{tt,cc} = a_0;
        MODELED_AMP{tt,cc}  = a_m;
        OBSERVED_PHS{tt,cc} = g_0;
        MODELED_PHS{tt,cc}  = g_m;
        V_cell{tt,cc} = V;
        Lat_cell{tt,cc} = T_o.Lat;
        Lon_cell{tt,cc} = T_o.Lon;
        
    end
end
error_matrix = D; 
error_values = D_cell; 
sta_locs     = {Lon_cell Lat_cell};

% Create Harm.txt file
fid=fopen('Harm.txt','w');  %<--this will be what this script produces!!
NoSta   = 0; % <--number of stations.
NoConst = length(const_cell); % <--number of constituents.
NoType  = length(type_cell); %<--different types of stations.

for CLASS = 1 : NoType
    NoSta = length(D_cell{CLASS,1}) + NoSta;
    NoStaPerClass(CLASS) = length(D_cell{CLASS,1});
end
fullmarkers={'s ','o ','h ','^ ','v ','x ','+ ','d '};
fprintf(fid,'%d %d\n',NoSta,NoConst) ; %<--change number of stations for type of stations
for ii = 1 : min(NoConst,8)
    marker = fullmarkers{ii};
    if ii == NoConst
        fprintf(fid,'%s\n',marker);
    else
        fprintf(fid,'%s',marker);  %<-- depending on length of const_cell, add more marker types.
    end
end
% for each type of station
StaCounter = 0; 
for CLASS = 1 : NoType
    % for each station
    for sta = 1 : NoStaPerClass(CLASS)
        StaCounter = StaCounter + 1; 
        fprintf(fid,'%s %f %f\n', num2str(StaCounter),Lon_cell{CLASS,1}(sta),Lat_cell{CLASS,1}(sta));
        % for each constitutent at each station (that exists). 
        for c = 1 : NoConst
            constituent = const_cell{c};
            if sta < length(OBSERVED_AMP{CLASS,c})
                %   const-name   obs-amp     mod-amp     obs-phs     mod-phs
                fprintf(fid,'%s %f %f %f %f\n',constituent,OBSERVED_AMP{CLASS,c}(sta),MODELED_AMP{CLASS,c}(sta),...
                    OBSERVED_PHS{CLASS,c}(sta),MODELED_PHS{CLASS,c}(sta));
            else
                % no data for this particular constituent.
                %   const-name   obs-amp     mod-amp     obs-phs     mod-phs
                fprintf(fid,'%s %f %f %f %f\n',constituent,NaN,NaN,NaN,NaN);
            end
        end
    end
end
fclose(fid);
end
