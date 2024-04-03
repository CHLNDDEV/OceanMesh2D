function write2dm(n, e, x, y, z, ikle, name)
    % Check if ikle is zero-based; if not, shift it to be zero-based
    if min(ikle(:)) == 1
        ikle = ikle - 1;
    end

    % Open file
    fid = fopen(name, 'w');
    
    % Write header
    fprintf(fid, 'MESH2D\n');
    
    % Prepare elements and nodes as cell arrays of strings
    elementLines = cell(e, 1);
    for i = 1:e
        elementLines{i} = sprintf('E3T %d %d %d %d 1\n', i, ikle(i, 1) + 1, ikle(i, 2) + 1, ikle(i, 3) + 1);
    end
    
    nodeLines = cell(n, 1);
    for i = 1:n
        nodeLines{i} = sprintf('ND %d %.6f %.6f %.6f\n', i, x(i), y(i), z(i));
    end
    
    % Combine and write all lines to file at once
    allLines = [elementLines; nodeLines];
    fprintf(fid, '%s', allLines{:});
    
    % Close file
    fclose(fid);
end
