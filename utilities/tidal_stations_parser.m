function obj = tidal_stations_parser(obj,sta_databases,type_cell)
% obj = tidal_stations_parser(obj,ts,te,sta_databases,type_cell)
% Input a msh class obj, read a string of sta_databases and puts the 
% locations and names into the f15 struct of the msh obj.
% type_cell is a cell of vectors where values 1, 2, and/or 3 correspond to 
% elevation, velocity, and met station outputs, respectively
% 
% The sta_database is a string corresponding to a url link that will get 
% stations from that database that fall within your mesh.  
% Currently just the XML of CO-OPS/NOS/NOAA stations are supported (set
% sta_database = 'CO-OPS').
%                                                                       
% Created by William Pringle March 16 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Getting the boundary of the mesh
%B  = boundary(obj.p(:,1),obj.p(:,2)); % can be slow 
B  = convhull(obj.p(:,1),obj.p(:,2));

nn = 0;
%% Loop over and read the station databases info 
for sta_database = sta_databases
    nn = nn + 1;
    type = type_cell{nn};
    if any(type == 1)
        % elevation
        obj.f15.nstae = 0; obj.f15.elvstaloc = []; obj.f15.elvstaname = [];
    end
    if any(type == 2)
        % velocity
        obj.f15.nstav = 0; obj.f15.velstaloc = []; obj.f15.velstaname = [];
    end
    if any(type == 3)
        % met
        obj.f15.nstam = 0; obj.f15.metstaloc = []; obj.f15.metstaname = [];
    end  
    if contains(sta_database,'csv')
        T = readtable(sta_database);
        Sta_num = height(T);
        Sta_lon = T.Lon; Sta_lat = T.Lat;
        Sta_type = type*ones(Sta_num,1);
        Sta_name = T.Name;
        Sta_ID   = strings(Sta_num,1);
    elseif strcmp(sta_database,'CO-OPS')
        xml = 'https://opendap.co-ops.nos.noaa.gov/stations/stationsXML.jsp';
        DOMnode = xmlread(xml);
        S = xml2struct(DOMnode);
        Sta_num = length(S.stations.station);
        Sta_lon = zeros(Sta_num,1); Sta_lat = zeros(Sta_num,1);
        Sta_type = zeros(Sta_num,3);
        Sta_name = strings(Sta_num,1); Sta_ID = strings(Sta_num,1);

        for s = 1:Sta_num
            Sta_lon(s) = str2double(S.stations.station{s}.metadata.location.long.Text);
            Sta_lat(s) = str2double(S.stations.station{s}.metadata.location.lat.Text);
            Sta_name{s} = S.stations.station{s}.Attributes.name;
            Sta_ID{s}   = S.stations.station{s}.Attributes.ID;
            if any(contains(fields(S.stations.station{s}),'parameter'))
                for i = 1:length(S.stations.station{s}.parameter)
                    % I don't think velocity info exists for these stations
                    if length(S.stations.station{s}.parameter) > 1
                        name = S.stations.station{s}.parameter{i}.Attributes.name;
                    else
                        name = S.stations.station{s}.parameter.Attributes.name;
                    end
                    if strcmp(name,'Water Level')
                        Sta_type(s,1) = 1;  
                    elseif strcmp(name,'Winds')
                        Sta_type(s,3) = 1;  
                    end
                end
            end
        end
    elseif strcmp(sta_database,'NDBC')
        Sta_num = 0; Sta_lon = 0; Sta_lat = 0; 
        Sta_type = zeros(1,3); Sta_name = strings; Sta_ID = strings;
        % Historical data
        HIST = urlread('http://www.ndbc.noaa.gov/historical_data.shtml'); 
        if any(type == 1)
            DART = urlread('http://www.ndbc.noaa.gov/dart_stations.php');
            C = textscan(DART,'%s','Delimiter','\t');
            C = C{1}; C(1:2) = [];
            for s = 1:length(C)
                Sta_ID(Sta_num + s) = C{s}(1:5);
                Sta_name(Sta_num + s) = strtrim(C{s}(14:64));
                Sta_lat(Sta_num + s) = str2double(C{s}(65:73));
                Sta_lon(Sta_num + s) = str2double(C{s}(76:85));
                Sta_type(Sta_num + s,1) = 1;
            end
            Sta_num = Sta_num + length(C);
            
            % water level
            ii1 = strfind(HIST,'<b>Water level');
            wl = HIST(ii1:end);
            ii = strfind(wl,'<li>');
            for s = 1:length(ii)
                Sta_ID(Sta_num + s) = wl(ii(s)+4:ii(s)+8);
                Sta_type(Sta_num + s,1) = 1; Sta_name(Sta_num + s) = '';
                Sta_lon(Sta_num + s) = 0; Sta_lat(Sta_num + s) = 0; 
            end   
            J = Sta_num + 1:Sta_num + length(ii);
            [Sta_lon(J),Sta_lat(J),Sta_name(J)] = read_NDBC_loc_name(Sta_ID(J));
            Sta_num = Sta_num + length(ii);
        end
        % Next first read NDBC stations
        if any(type > 1) 
            NDBC = urlread('http://www.ndbc.noaa.gov/stndesc.shtml');
            ii1 = strfind(NDBC,'<pre>'); ii2 = strfind(NDBC,'</pre>');
            NDBC = NDBC(ii1:ii2);
            C = textscan(NDBC,'%s','Delimiter','\t');
            C = C{1}; C(1:5) = []; C(end) = [];
            for s = 1:length(C)
                NDBC = textscan(C{s},'%d %s %f %f %f %s %f %f %f %s %d %d %s');
                NDBC_ID(s) = string(NDBC{1}(1));
                NDBC_name(s) = string(NDBC{2});
                N = NDBC{6}; E = NDBC{10};
                NDBC_lat(s) = NDBC{3} + NDBC{4}/60 + NDBC{5}/3600;
                if strcmp(N{1}(1),'S')
                    NDBC_lat(s) = -NDBC_lat(s);
                end
                NDBC_lon(s) = double(NDBC{7} + NDBC{8}/60 + NDBC{9}/3600);
                if strcmp(E{1}(1),'W')
                    NDBC_lon(s) = -NDBC_lon(s);
                end
            end
            if any(type == 2)
                % Ocean current 
                ii2 = strfind(HIST,'<b>Ocean Current Data');
                ii1 = strfind(HIST,'<b>Oceanographic Data');
                wl = HIST(ii2:ii1-35);
                ii = strfind(wl,'<li>');
                for s = 1:length(ii)
                    Sta_ID(Sta_num + s) = wl(ii(s)+4:ii(s)+8);
                    Sta_type(Sta_num + s,2) = 1; Sta_name(Sta_num + s) = '';
                    Sta_lon(Sta_num + s) = 0; Sta_lat(Sta_num + s) = 0; 
                end 
                J = Sta_num + 1:Sta_num + length(ii);
                [~,K,M] = intersect(Sta_ID(J),NDBC_ID);
                Sta_lon(J(K)) = NDBC_lon(M); Sta_lat(J(K)) = NDBC_lat(M);
                Sta_name(J(K)) = NDBC_name(M);
                [~,K] = setdiff(Sta_ID(J),NDBC_ID);
                [Sta_lon(J(K)),Sta_lat(J(K)),Sta_name(J(K))] = read_NDBC_loc_name(Sta_ID(J(K)));
                Sta_num = Sta_num + length(ii);
            end
            if any(type == 3)
                % meteorological, copy over NDBC info
                J = Sta_num + 1:Sta_num + length(NDBC_ID);
                Sta_ID(J) = NDBC_ID;
                Sta_lon(J) = NDBC_lon;
                Sta_lat(J) = NDBC_lat;
                Sta_name(J) = NDBC_name;
                Sta_type = [Sta_type; zeros(length(J),3)]; 
                Sta_type(J,3) = 1;
                Sta_num = Sta_num + length(NDBC_ID);
            end
        end
        Sta_lon = Sta_lon';  Sta_lat = Sta_lat'; 
        Sta_ID = Sta_ID'; Sta_name = Sta_name';
        
    else
        error([char(sta_database) ' not yet supported']) 
    end

    %% Delete stuff outside the boundary of the mesh
    if ~contains(sta_database,'csv')
        in = inpoly([Sta_lon Sta_lat],obj.p(B,:));
        Sta_lon(~in) = []; Sta_lat(~in) = []; Sta_type(~in,:) = [];
        Sta_name(~in) = []; Sta_ID(~in) = [];
    end
    
    %% Now put into f15 struct
    if any(type == 1)
        % elevation
        obj.f15.nstae = obj.f15.nstae + numel(find(Sta_type(:,1) == 1));
        obj.f15.elvstaloc = [obj.f15.elvstaloc;
                             [Sta_lon(Sta_type(:,1) == 1) ...
                             Sta_lat(Sta_type(:,1) == 1)]];
        obj.f15.elvstaname = [obj.f15.elvstaname;
                             strcat(Sta_name(Sta_type(:,1) == 1),' ID:',...
                                    Sta_ID(Sta_type(:,1) == 1))];
    end
    if any(type == 2)
        % velocity
        obj.f15.nstav = obj.f15.nstav + numel(find(Sta_type(:,2) == 1));
        obj.f15.velstaloc = [obj.f15.velstaloc;
                             [Sta_lon(Sta_type(:,2) == 1) ...
                             Sta_lat(Sta_type(:,2) == 1)]];
        obj.f15.velstaname = [obj.f15.velstaname;
                             strcat(Sta_name(Sta_type(:,2) == 1),' ID:',...
                                    Sta_ID(Sta_type(:,2) == 1))];
    end
    if any(type == 3)
        % met
        obj.f15.nstam = obj.f15.nstam + numel(find(Sta_type(:,3) == 1));
        obj.f15.metstaloc = [obj.f15.metstaloc;
                             [Sta_lon(Sta_type(:,3) == 1) ...
                             Sta_lat(Sta_type(:,3) == 1)]];
        obj.f15.metstaname = [obj.f15.metstaname;
                             strcat(Sta_name(Sta_type(:,3) == 1),' ID:',...
                                    Sta_ID(Sta_type(:,3) == 1))];
    end    
end

end

function [lon,lat,name] = read_NDBC_loc_name(Sta_ID)
    Prefix = 'http://www.ndbc.noaa.gov/station_page.php?station=';
    lat = zeros(length(Sta_ID),1);  lon = zeros(length(Sta_ID),1);
    name = strings(length(Sta_ID),1);
    for s = 1:length(Sta_ID)
        NDBC = urlread([Prefix Sta_ID{s}]);
        ii = strfind(NDBC,'('); wl = NDBC(ii(1)+1:ii(1)+100);
        C = textscan(wl,'%f %s %f %s');
        N = C{2}; E = C{4};
        lat(s) = C{1};
        if strcmp(N{1}(1),'S')
            lat(s) = -lat(s);
        end
        lon(s) = C{3};
        if strcmp(E{1}(1),'W')
            lon(s) = -lon(s);
        end
        ii1 = strfind(wl,')'); ii2 = strfind(wl,'."'); 
        name(s) = wl(ii1+4:ii2-1);
    end
end