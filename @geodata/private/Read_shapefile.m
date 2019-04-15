function polygon_struct = Read_shapefile( finputname, polygon, bbox, ...
                                          h0, boubox, plot_on )
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
    % bbox straddles 180/-180 line
    loop = 2; minus = 0;
elseif bbox(1,1) > 180 && bbox(1,1) > 180
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
            % The shaperead is much faster if it is available
            % Read the structure
            try 
                S = shaperead(fname{1},'BoundingBox',bboxt);
                % Get rid of unwanted components;
                D = struct2cell(S);
                S = cell2struct(D(3:4,:)',{'X','Y'},2);
                disp('Read shapefile with shaperead')
                sr = 1;
            catch
                disp('Reading shapefile with m_shaperead')
                % This uses m_map (slower but free)
                S = m_shaperead(fname{1},reshape(bboxt',4,1));
                % Let's just keep the x-y data
                D = S.ncst;
                S = cell2struct(D','points',1);
                sr = 0;
            end
            if ~isempty(S)
                % Keep the following polygons
                SG = [SG; S];
            end
        end
    end
else
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
    tmpCC = []; nn = 0;
    for ii = 1:length(tmpC)
       % may have NaNs inside 
       isnan1 = find(isnan(tmpC{ii}(:,1)));
       if isempty(isnan1)
           isnan1 = length(tmpC{ii})+1; 
       elseif isnan1(end) ~= length(tmpC{ii})
           isnan1(end+1) = length(tmpC{ii})+1;
       end
       is = 1;
       for jj = 1:length(isnan1)
           nn = nn + 1;
           ie = isnan1(jj)-1;
           tmpCC{nn} = tmpC{ii}(is:ie,:);
           is = isnan1(jj)+1;
       end
    end
    if ~isempty(tmpCC)
        tmpC = tmpCC';
    end  
    tmpM =  cell2mat(tmpC); 
    if size(tmpM,2) == 3
       tmpM = tmpM(:,1:2); 
    end
end
% Get current polygon
% Check proportion of polygon that is within bbox
tmpIn = inpoly(tmpM,polygon_struct.outer, edges);
tmpInC = mat2cell(tmpIn,cellfun(@length,tmpC));

j = 0 ; k = 0 ; height = []; new_islandb = []; new_mainb = [];
for i = 1 : length(tmpC)
    if sr
        points = tmpC{i}(1:end-1,:) ;
        In     = tmpInC{i}(1:end-1) ;
    else
        points = tmpC{i}(1:end,:) ;
        if size(points,2) == 3
            height = points(:,3); 
            points = points(:,1:2); 
        else
            height = [];
        end
        In     = tmpInC{i}(1:end) ;
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
        end
    end
end
if k > 0
    polygon_struct.inner = [polygon_struct.inner; cell2mat(new_island')];
    polygon_struct.innerb = cell2mat(new_islandb');
end
if j > 0
    polygon_struct.mainland = [polygon_struct.mainland; cell2mat(new_main')];
    polygon_struct.mainlandb = cell2mat(new_mainb');
end
% Merge overlapping mainland and inner
if ~isempty(new_mainb) & k > 0
    polym = polygon_struct.mainland;
    mergei = false(k,1);
    for kk = 1:k
        polyi = new_island{kk};
        IA = find(ismembertol(polym,polyi,1e-5,'ByRows',true));
        if length(IA) > 2
           [x,y] = polybool('or',polym(:,1),...
                                        polym(:,2),polyi(:,1),polyi(:,2));
           polym = [x,y];
           mergei(kk) = 1;
        end
    end
    if ~isnan(polym(end,1)); polym(end+1,:) = NaN; end
    polygon_struct.mainland = polym;
    new_island(mergei) = [];
    polygon_struct.inner = cell2mat(new_island');
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
