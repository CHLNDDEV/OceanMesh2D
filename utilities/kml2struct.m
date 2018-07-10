function kmlStruct = kml2struct(kmlFile)
% kmlStruct = kml2struct(kmlFile)
%
% Import a .kml file as a vector array of shapefile structs, with Geometry, Name,
% Description, Lon, Lat, and BoundaryBox fields.  Structs may contain a mix
% of points, lines, and polygons.
%
% .kml files with folder structure will not be presented as such, but will
% appear as a single vector array of structs.
%
% 

[FID msg] = fopen(kmlFile,'rt');

if FID<0
    error(msg)
end

txt = fread(FID,'uint8=>char')';
fclose(FID);

expr = '<Placemark.+?>.+?</Placemark>';

objectStrings = regexp(txt,expr,'match');

Nos = length(objectStrings);

for ii = 1:Nos
    % Find Object Name Field
    bucket = regexp(objectStrings{ii},'<name.*?>.+?</name>','match');
    if isempty(bucket)
        name = 'undefined';
    else
        % Clip off flags
        name = regexprep(bucket{1},'<name.*?>\s*','');
        name = regexprep(name,'\s*</name>','');
        if strcmp(name(1:2),'<!')
            name = regexprep(name,'<![CDATA[\s*','');
            name = regexprep(name,']]>\s*','');
        end
    end
    
    % Find Object Description Field
    bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
    if isempty(bucket)
        desc = '';
    else
        % Clip off flags
        desc = regexprep(bucket{1},'<description.*?>\s*','');
        desc = regexprep(desc,'\s*</description>','');
        
        % Clip off other flags and rearrange into rows
        if strcmp(desc(1:2),'<!')
            desc = regexprep(desc,'<![CDATA[\s*','');
            desc = regexprep(desc,']]>\s*','');
            % Splitting into rows
            I = find(desc == '<');
            if ~isempty(I)
                desc_c = cell(length(I)+1,1);
                desc_c{1} = desc(1:I(1)-1);
                for i = 2:length(I)
                    desc_c{i} = desc(I(i-1)+4:I(i)-1);
                end
                desc_c{end} = desc(I(end)+4:end);
            end
        end
    end
        
    % Find style
    bucket = regexp(objectStrings{ii},'<styleUrl.*?>.+?</styleUrl>','match');
    if isempty(bucket)
        StyleUrl = 'undefined';
    else
        % Clip off flags
        StyleUrl = regexprep(bucket{1},'<styleUrl.*?>\s*','');
        StyleUrl = regexprep(StyleUrl,'\s*</styleUrl>','');
    end
    
    geom = 0;
    % Identify Object Type
    if ~isempty(regexp(objectStrings{ii},'<Point', 'once'))
        geom = 1;
    elseif ~isempty(regexp(objectStrings{ii},'<LineString', 'once'))
        geom = 2;
    elseif ~isempty(regexp(objectStrings{ii},'<Polygon', 'once'))
        geom = 3;
    end
    
    switch geom
        case 1
            geometry = 'Point';
        case 2
            geometry = 'Line';
        case 3
            geometry = 'Polygon';
        otherwise
            geometry = '';
    end
    
    % Find Coordinate Field
    bucket = regexp(objectStrings{ii},'<coordinates.*?>.+?</coordinates>','match');
    % Clip off flags
    coordStr = regexprep(bucket{1},'<coordinates.*?>(\s+)*','');
    coordStr = regexprep(coordStr,'(\s+)*</coordinates>','');
    % Split coordinate string by commas or white spaces, and convert string
    % to doubles
    coordMat = str2double(regexp(coordStr,'[,\s]+','split'));
    % Rearrange coordinates to form an x-by-3 matrix

    [m,n] = size(coordMat);
    if m*n >= 3
        coordMat = reshape(coordMat,3,m*n/3)';
    end
    
    % define polygon in clockwise direction, and terminate
    [Lat, Lon] = poly2ccw(coordMat(:,2),coordMat(:,1));
    if geom==3
        Lon = [Lon;NaN];
        Lat = [Lat;NaN];
    end
    
    % Create structure
    kmlStruct(ii).Geometry = geometry;
    kmlStruct(ii).Name = name;
    if exist('desc_c','var')
        kmlStruct(ii).Description = desc_c;
    else
        kmlStruct(ii).Description = desc;
    end
    kmlStruct(ii).Lon = Lon;
    kmlStruct(ii).Lat = Lat;
    kmlStruct(ii).BoundingBox = [[min(Lon) min(Lat);max(Lon) max(Lat)]];
    kmlStruct(ii).StyleUrl = StyleUrl;
end