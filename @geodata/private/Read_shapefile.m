function polygon_struct = Read_shapefile( finputname, polygon, bbox, ...
    h0, plot_on )
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
% OUTPUTS: 
% polygon_struct    : a structure containing the vector features identified as
%              either islands, mainland, or outer.
% Written by William Pringle and Keith Roberts, CHL,UND, 2017
%% Loop over all the filenames and get the shapefile within bbox
SG = [];
if (size(finputname,1)~=0)
    for fname = finputname
        % The shaperead is much faster if it is available
        if exist('shaperead','file')
            disp('Reading coastline with shaperead')
            % Read the structure
            S = shaperead(fname{1},'BoundingBox',bbox');
            % Get rid of unwanted components;
            D = struct2cell(S);
            S = cell2struct(D(3:4,:)',{'X','Y'},2);
        else
            disp('Reading coastline with m_shaperead')
            % This uses m_map (slower but free)
            S = m_shaperead(fname{1},bbox(:)');
            % Let's just keep the x-y data
            D = S.ncst;
            S = cell2struct(D','points',1);
        end
        if ~isempty(S) 
            % Keep the following polygons
            SG = [SG; S];
        end
    end
else
    count = 1;
    j=1;
    for i = 1 : length(polygon)
        if(isnan(polygon(i,1))==1)
            count = count + 1; j=1;
            continue 
        end
        SG(count).points(j,:) = polygon(i,:);
        %SG(count).X(j) = polygon(i,1);
        %SG(count).Y(j) = polygon(i,2);
        j=j+1; 
    end
end
% If we don't have an outer polygon already then make it by bbox
polygon_struct.outer = [bbox(1,1) bbox(2,1);
    bbox(1,1) bbox(2,2);
    bbox(1,2) bbox(2,2);
    bbox(1,2) bbox(2,1);
    bbox(1,1) bbox(2,1)];
% Densify the outer polygon (fills gaps larger than half min edgelength).
%[latout,lonout] = interpm(polygon_struct.outer(:,2),...
%    polygon_struct.outer(:,1),h0/2);
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
edges = Get_poly_edges( [polygon_struct.outer; NaN NaN] );

for s = 1:length(SG)
    if exist('shaperead','file')
        points = [SG(s).X' SG(s).Y'];
        points(isnan(points(:,1)),:) = [];
    else
        points = SG(s).points;
    end
    % Get current polygon
    % Check proportion of polygon that is within bbox
    In = inpoly(points,[polygon_struct.outer; NaN NaN], edges);

    % lets calculate the area of the
    % feature using the shoelace algorithm and decided whether to keep or
    % not based on area.
    if all(points(1,:) == points(end,:))
        area = shoelace(points(:,1),points(:,2));
    else
        area = 999; % not a polygon
    end
    if (area < 4*h0^2) % too small, then don't consider it.
        continue;
    elseif length(find(In == 1)) == length(points)
        % Wholey inside box, set as island (with NaN delimiter)
        new_island = [points; NaN NaN];
        polygon_struct.inner = [polygon_struct.inner; new_island];
    else
        % Partially inside box, set as mainland
        new_main = [points; NaN NaN];
        polygon_struct.mainland = [polygon_struct.mainland; new_main];
    end
end

polygon_struct.outer = [polygon_struct.outer; NaN NaN; polygon_struct.mainland];

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
