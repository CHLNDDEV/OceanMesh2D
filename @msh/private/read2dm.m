function [n, e, x, y, z, ikle] = read2dm(two_dm_file)
    % Validate the file
    if ~isfile(two_dm_file)
        error('File does not exist.');
    end

    % Open file
    fid = fopen(two_dm_file, 'r');

    % Skip header
    fgetl(fid);

    % Initialize
    n = 0;
    e = 0;
    x = [];
    y = [];
    z = [];
    ikle = [];

    % Read in chunks
    chunkSize = 1e5;  % Adjust based on your system's memory
    while ~feof(fid)
        chunk = textscan(fid, '%s', chunkSize, 'Delimiter', '\n');
        chunk = chunk{1};

        % Pre-allocate for speed
        n_tmp = 0;
        e_tmp = 0;
        x_tmp = zeros(length(chunk), 1);
        y_tmp = zeros(length(chunk), 1);
        z_tmp = zeros(length(chunk), 1);
        ikle_tmp = zeros(length(chunk), 3);

        for i = 1:length(chunk)
            line = chunk{i};
            if startsWith(line, 'E3T')
                e_tmp = e_tmp + 1;
                tokens = sscanf(line, 'E3T %*d %d %d %d');
                ikle_tmp(e_tmp, :) = tokens' - 1;
            elseif startsWith(line, 'ND')
                n_tmp = n_tmp + 1;
                tokens = sscanf(line, 'ND %*d %f %f %f');
                x_tmp(n_tmp) = tokens(1);
                y_tmp(n_tmp) = tokens(2);
                z_tmp(n_tmp) = tokens(3);
            end
        end

        % Trim and append
        x_tmp = x_tmp(1:n_tmp);
        y_tmp = y_tmp(1:n_tmp);
        z_tmp = z_tmp(1:n_tmp);
        ikle_tmp = ikle_tmp(1:e_tmp, :);

        x = [x; x_tmp];
        y = [y; y_tmp];
        z = [z; z_tmp];
        ikle = [ikle; ikle_tmp];

        n = n + n_tmp;
        e = e + e_tmp;
    end

    fclose(fid);
end
