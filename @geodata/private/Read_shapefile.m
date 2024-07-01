function polygon_struct = Read_shapefile( finputname, polygon, bbox, ...
    h0, boubox, plot_on, shapefile_3d)
% Read_shapefile: Reads a shapefile or a NaN-delimited vector
% containing polygons and/or segments in the the desired region
% of interest. Classifies the vector data as either a
% mainland, island, or outer boundary. The program will automatically
% trim islands that are have an area smaller than 4*h0^2 and will in gaps
% in the vector that are larger than h0/2. This is necessary for the
% use of boundary rejection method in dpoly.

% INPUTS:
% finputname : file name(s) of the shapefile listed in a cell array
% polygon    : a NaN-delimited vector of polygons and/or segments.
% bbox       : the bounding box that we want to extract out
% h0         : minimum edge length in region of interest.
% plot_on    : plot the final polygon or not (=1 or =0)
%
% OUTPUTS:
% polygon_struct    : a structure containing the vector features identified as
%              either islands, mainland, or outer.
% Written by William Pringle and Keith Roberts, CHL,UND, 2017
% Edits by Keith Roberts, July 2018.
%% Loop over all the filenames and get the shapefile within bbox
SG = [];
if bbox(1,2) > 180 && bbox(1,1) < 180
    % bbox straddles 180/-180 meridian
    loop = 2; minus = 0;
elseif all(bbox(1,:) > 180)
    % beyond 180 in 0 to 360 format
    loop = 1; minus = 1;
else
    loop = 1; minus = 0;
end
if (size(finputname,1)~=0)
    for fname = finputname
        for nn = 1:loop
            bboxt = bbox';
            if loop == 2
                if nn == 1
                    bboxt(2,1) = 180;
                else
                    bboxt(1,1) = -180; bboxt(2,1) = bboxt(2,1) - 360;
                end
            end
            if minus
                bboxt(:,1) = bboxt(:,1) - 360;
            end
            % Read the structure
            try
                % The shaperead is faster if it is available
                S = shaperead(fname{1},'BoundingBox',bboxt);
                % Get rid of unwanted components;
                D = struct2cell(S);
                S = cell2struct(D(3:4,:)',{'X','Y'},2);
                disp('Read shapefile with shaperead')
                sr = 1;
            catch
                % If only m_shaperead is available or if some error occured
                % with shaperead (e.g., 3D shapefile, m_shaperead may work)
                disp('Reading shapefile with m_shaperead')
                % This uses m_map (slower but free)
                S = m_shaperead(fname{1},reshape(bboxt',4,1));
                % Let's just keep the x-y data
                D = S.ncst;
                if isfield(S,'dbf')  || isfield(S,'dbfdata')
                    code = S.dbfdata(:,1);
                    S = cell2struct([D code]',{'points' 'type'},1);
                else
                    S = cell2struct(D','points',1);
                end
                sr = 0;
            end
            if ~isempty(S)
                % Keep the following polygons
                SG = [SG; S];
            end
        end
    end
else
    sr = 1 ;
    % convert NaN-delimited vector to struct
    count = 1; j=1;
    for i = 1 : length(polygon)
        % the end of the segment
        
        if(isnan(polygon(i,1))==1)
            % put NaN at end
            SG(count,1).X(:,j) =NaN;
            SG(count,1).Y(:,j) =NaN;
            % reset
            j = 1 ; count = count + 1;
            
            continue
        else
            % keep going
            SG(count,1).X(:,j) = polygon(i,1);
            SG(count,1).Y(:,j) = polygon(i,2);
            j=j+1;
        end
    end
end
% If we don't have an outer polygon already then make it by bbox
polygon_struct.outer = boubox;
% Densify the outer polygon (fills gaps larger than half min edgelength).
[latout,lonout] = my_interpm(polygon_struct.outer(:,2),...
    polygon_struct.outer(:,1),h0/2);
polygon_struct.outer = [];
polygon_struct.outer(:,1) = lonout;
polygon_struct.outer(:,2) = latout;

%% Find whether the polygons are wholey inside the bbox..
%% Set as islands or mainland
disp('Partitioning the boundary into islands, mainland, ocean')
polygon_struct.inner = [];
polygon_struct.mainland = [];
polygon_struct.innerb = [];
polygon_struct.mainlandb = [];
polygon_struct.innerb_type = [];
polygon_struct.mainlandb_type = [];
edges = Get_poly_edges( polygon_struct.outer );

if isempty(SG); return; end

if sr
    tmpM = [[SG.X]',[SG.Y]'] ; % MAT
    if bbox(1,2) > 180
        tmpM(tmpM(:,1) < 0,1) =  tmpM(tmpM(:,1) < 0,1) + 360;
    end
    for i = 1 : length(SG)
        dims(i) = length(SG(i).X) ;
    end
    tmpC = mat2cell(tmpM,dims); % TO CELL
else
    tmpC =  struct2cell(SG)';
    
    % Remove empty cells from shapefile
    if any(cellfun(@isempty,tmpC(:,1)))
        warning('Empty objects in shapefile have been removed')
        tmpC = tmpC(~cellfun(@isempty,tmpC(:,1)),:);
    end

    tmpCC = []; nn = 0;
    for ii = 1:size(tmpC,1)
        % may have NaNs inside
        isnan1 = find(isnan(tmpC{ii,1}(:,1)));
        if isempty(isnan1)
            isnan1 = length(tmpC{ii,1})+1;
        elseif isnan1(end) ~= length(tmpC{ii,1})
            isnan1(end+1) = length(tmpC{ii,1})+1;
        end
        is = 1;
        for jj = 1:length(isnan1)
            nn = nn + 1;
            ie = isnan1(jj)-1;
            tmpCC{nn,1} = tmpC{ii,1}(is:ie,:);
            tmpCC{nn,2} = tmpC{ii,2};
            is = isnan1(jj)+1;
        end
    end
    if ~isempty(tmpCC)
        tmpC = tmpCC;
    end
    tmpM =  cell2mat(tmpC(:,1));
    if size(tmpM,2) == 3
        tmpM = tmpM(:,1:2);
    end
end
% Get current polygon
% Check proportion of polygon that is within bbox
tmpIn = inpoly(tmpM,polygon_struct.outer, edges);
tmpInC = mat2cell(tmpIn,cellfun(@length,tmpC(:,1)));

j = 0 ; k = 0 ; height = []; new_islandb = []; new_mainb = [];
new_islandb_type = []; new_mainb_type = [];
for i = 1 : size(tmpC,1)
    if sr
        % using shaperead
        points = tmpC{i,1}(1:end-1,:) ;
        In     = tmpInC{i,1}(1:end-1) ;
    else
        % using m_shaperead
        points = tmpC{i,1}(1:end,:) ;
        if size(points,2) == 3
            points = points(:,1:2);
        end
        if shapefile_3d
            % if 3-D shapefile
            height = points(:,3);
            points = points(:,1:2);
            type   = tmpC{i,2};
            if strcmp(type,'BA040')
                type = 'ocean';
            elseif strcmp(type,'BH080')
                type = 'lake';
            elseif strcmp(type,'BH140')
                type = 'river';
            end
        else
            height = [];
        end
        In     = tmpInC{i,1}(1:end) ;
    end
    if bbox(1,2) > 180
        lond = abs(diff(points(:,1)));
        if any(lond > 350)
            points(points(:,1) > 180,1) = 0;
        end
    end
    % lets calculate the area of the
    % feature using the shoelace algorithm and decided whether to keep or
    % not based on area.
    if all(points(1,:) == points(end,:))
        area = shoelace(points(:,1),points(:,2));
    else
        area = 999; % not a polygon
    end
    if length(find(In == 1)) == length(points)
        % Wholey inside box
        if area < 4*h0^2 % too small, then don't consider it.
            continue;
        end
        % Set as island (with NaN delimiter)
        k = k + 1 ;
        new_island{k} = [points; NaN NaN];
        if ~isempty(height)
            new_islandb{k} = [points height; NaN NaN NaN];
            new_islandb_type{k} = type;
        end
    else
        %Partially inside box
        if area < 100*h0^2 % too small, then don't consider it.
            continue;
        end
        % Set as mainland
        j = j + 1 ;
        new_main{j} = [points; NaN NaN];
        if ~isempty(height)
            new_mainb{j} = [points height; NaN NaN NaN];
            new_mainb_type{j} = type;
        end
    end
end
if k > 0
    polygon_struct.inner = [polygon_struct.inner; cell2mat(new_island')];
    polygon_struct.innerb = cell2mat(new_islandb');
    polygon_struct.innerb_type = new_islandb_type;
end
if j > 0
    polygon_struct.mainland = [polygon_struct.mainland; cell2mat(new_main')];
    polygon_struct.mainlandb = cell2mat(new_mainb');
    polygon_struct.mainlandb_type = new_mainb_type;
end
% Merge overlapping mainland and inner
if exist('polybool','file') || exist('polyshape','file')
    if ~isempty(new_mainb) && k > 0
        polym = polygon_struct.mainland;
        mergei = false(k,1);
        for kk = 1:k
            polyi = new_island{kk};
            IA = find(ismembertol(polym,polyi,1e-5,'ByRows',true));
            if length(IA) > 2
                if exist('polyshape','file')
                    polyout = union(polyshape(polym),polyshape(polyi));
                    polym = polyout.Vertices;
                else
                    [x,y] = polybool('union',polym(:,1),polym(:,2),...
                        polyi(:,1),polyi(:,2));
                    polym = [x,y];
                end
                mergei(kk) = 1;
            end
        end
        if ~isnan(polym(end,1)); polym(end+1,:) = NaN; end
        polygon_struct.mainland = polym;
        new_island(mergei) = [];
        polygon_struct.inner = cell2mat(new_island');
    end
else
    warning(['no polyshape or polybool available to merge possible ' ...
        'overlapping of mainland and inner'])
end
% Remove parts of inner and mainland overlapping with outer
polygon_struct.outer = [polygon_struct.outer; polygon_struct.mainland];
%% Plot the map
if plot_on >= 1 && ~isempty(polygon_struct)
    figure(1);
    hold on
    plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2))
    if ~isempty(polygon_struct.inner)
        plot(polygon_struct.inner(:,1),polygon_struct.inner(:,2))
    end
    if ~isempty(polygon_struct.mainland)
        plot(polygon_struct.mainland(:,1),polygon_struct.mainland(:,2))
    end
end
%EOF
end
