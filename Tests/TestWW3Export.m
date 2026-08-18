function TestWW3Export
% Regression test for ADCIRC -> WW3 export with sparse open boundaries.

fort14 = [tempname '.14'];
ww3base = tempname;
ww3file = [ww3base '.ww3'];
ww3emptybase = tempname;
ww3emptyfile = [ww3emptybase '.ww3'];
cleanupObj = onCleanup(@() cleanup_temp_files({fort14, ww3file, ww3emptyfile})); %#ok<NASGU>

write_temp_fort14(fort14);

m = msh('fname',fort14);
write(m,ww3base,'ww3');

if exist(ww3file,'file') ~= 2
    error('WW3 export regression failed: output file was not created.');
end

lines = read_text_lines(ww3file);
[element_count, element_lines] = get_element_block(lines);
expected_boundary_nodes = [1; 3; 4];
expected_triangle_count = size(m.t,1);
expected_element_count = expected_triangle_count + length(expected_boundary_nodes);

if element_count ~= expected_element_count
    error('Incorrect WW3 element count for sparse boundary export. Got %d, expected %d.', ...
        element_count, expected_element_count);
end

boundary_nodes = parse_boundary_nodes(element_lines, length(expected_boundary_nodes));
if ~isequal(boundary_nodes(:), expected_boundary_nodes)
    error('Incorrect WW3 boundary nodes. Got [%s], expected [%s].', ...
        num2str(boundary_nodes(:)'), num2str(expected_boundary_nodes'));
end

if any(boundary_nodes <= 0)
    error('WW3 export wrote invalid padded boundary node IDs.');
end

m.op = [];
write(m,ww3emptybase,'ww3');

if exist(ww3emptyfile,'file') ~= 2
    error('WW3 empty-boundary regression failed: output file was not created.');
end

lines_empty = read_text_lines(ww3emptyfile);
[element_count_empty, element_lines_empty] = get_element_block(lines_empty);
if element_count_empty ~= expected_triangle_count
    error('Incorrect WW3 element count for empty-boundary export. Got %d, expected %d.', ...
        element_count_empty, expected_triangle_count);
end

if length(element_lines_empty) ~= expected_triangle_count
    error('Unexpected number of WW3 element lines for empty-boundary export.');
end

fprintf('Passed: WW3 sparse-boundary export regression\n');

end

function write_temp_fort14(fort14)

fid = fopen(fort14,'w');
if fid < 0
    error('Could not create temporary fort.14 file for WW3 regression test.');
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'WW3 export regression mesh\n');
fprintf(fid,'2 4\n');
fprintf(fid,'1 0 0 -1\n');
fprintf(fid,'2 1 0 -2\n');
fprintf(fid,'3 1 1 -3\n');
fprintf(fid,'4 0 1 -4\n');
fprintf(fid,'1 3 1 2 3\n');
fprintf(fid,'2 3 1 3 4\n');
fprintf(fid,'2\n');
fprintf(fid,'3\n');
fprintf(fid,'2\n');
fprintf(fid,'1\n');
fprintf(fid,'3\n');
fprintf(fid,'1\n');
fprintf(fid,'4\n');
fprintf(fid,'0\n');
fprintf(fid,'0\n');

end

function lines = read_text_lines(fname)

raw = fileread(fname);
lines = regexp(raw,'\r\n|\n|\r','split');
if ~isempty(lines) && isempty(lines{end})
    lines(end) = [];
end

end

function [element_count, element_lines] = get_element_block(lines)

idx_start = find(strcmp(lines,'$Elements'),1,'first');
idx_end = find(strcmp(lines,'$EndElements'),1,'first');
if isempty(idx_start) || isempty(idx_end) || idx_end <= idx_start + 1
    error('Could not locate a valid WW3 $Elements block in the output file.');
end

element_count = sscanf(lines{idx_start + 1},'%d',1);
element_lines = lines(idx_start + 2:idx_end - 1);

if length(element_lines) ~= element_count
    error('WW3 element block length does not match the declared element count.');
end

end

function boundary_nodes = parse_boundary_nodes(element_lines, nb_boundary_nodes)

boundary_nodes = zeros(nb_boundary_nodes,1);
for i = 1:nb_boundary_nodes
    values = sscanf(element_lines{i},'%d');
    if length(values) < 6
        error('Malformed WW3 boundary element line: %s', element_lines{i});
    end
    if values(2) ~= 15
        error('Expected WW3 point element type 15 for boundary element %d.', i);
    end
    boundary_nodes(i) = values(end);
end

end

function cleanup_temp_files(files)

for i = 1:length(files)
    if exist(files{i},'file') == 2
        delete(files{i});
    end
end

end
