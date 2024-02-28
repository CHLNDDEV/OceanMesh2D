function polygons = nandelim_to_cell(data)
%NANDELIM_TO_CELL Convert NaN-delimited data into a cell array of polygons.
%   polygons = NANDELIM_TO_CELL(data) takes an Nx2 matrix 'data', where
%   each polygon is separated by a row with NaN in the first column, and
%   converts it into a cell array. Each cell contains the vertices of a
%   polygon defined by consecutive rows up to the next NaN row.

    % Identify the indices of rows that have NaN in the first column
    nanRows = find(isnan(data(:, 1)));

    % Preallocate the cell array for efficiency
    polygons = cell(length(nanRows) + 1, 1);

    % Set the starting index for the first polygon
    startIndex = 1;

    % Iterate over each NaN row to extract polygons
    for i = 1:length(nanRows)
        % Extract the polygon data up to the row before the NaN
        polygons{i} = data(startIndex:nanRows(i) - 1, :);

        % Update the start index for the next polygon
        startIndex = nanRows(i) + 1;
    end

    % Handle the last polygon after the final NaN row
    polygons{end} = data(startIndex:end, :);
end
