function area = shoelace(x, y)
    % Computes the area of a polygon using the Shoelace formula.
    % Handles NaNs by computing areas of individual sub-polygons separately.

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);

    % Identify NaN indices to split polygons
    nanIndices = isnan(x) | isnan(y);
    
    % Find segment start and end indices
    segmentStart = find([true; nanIndices(1:end-1)] & ~nanIndices);
    segmentEnd = find(~nanIndices & [nanIndices(2:end); true]);

    % Initialize total area
    area = 0;

    % Process each sub-polygon separately
    for i = 1:length(segmentStart)
        xi = x(segmentStart(i):segmentEnd(i));
        yi = y(segmentStart(i):segmentEnd(i));

        % Ensure the polygon is closed (first and last points should match)
        if ~isequal(xi(1), xi(end)) || ~isequal(yi(1), yi(end))
            xi = [xi; xi(1)];
            yi = [yi; yi(1)];
        end

        % Apply Shoelace formula
        Ai = 0.5 * abs(sum(xi(1:end-1) .* yi(2:end) - xi(2:end) .* yi(1:end-1)));
        area = area + Ai; % Accumulate area from all sub-polygons
    end
end
