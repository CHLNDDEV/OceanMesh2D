function fh = enforce_min_ef(fh)
% Ensures that the minimum mesh sizing value is always used in the case 
% of many edgefxs (e.g., multiscale meshing technique).
%
% Loop through edgefunctions to check if they overlap. 
% If they overlap, find the points that lie within the 
% upper edgefunction and compare them to the ones below.
%
% By Coleman Blakely UND, 2020.
%   Modifications by Keith Roberts USP, 2020. 

for kk = length(fh):-1:2
    fh_test = fh{kk};
    % make the points to use with inpoly to find the overlapping points
    [test_pts, ~] = edgefx2pts(fh_test);
    test_vals = fh_test.F.Values;
    x_test = fh_test.F.GridVectors{1};
    y_test = fh_test.F.GridVectors{2};
    [X_test, Y_test] = ndgrid(x_test,y_test);
    % Loop through all boxes that may be below the test box
    for ii = kk-1:-1:1
        fh_lower = fh{ii};
        % find the test points that are within the lower boundary
        id_overlap = inpoly(test_pts,fh_lower.boubox(1:end-1,:));
        % check to see if the test box overlaps the lower box. If not, skip
        % this box
        if sum(id_overlap)==0
            continue
        end
        overlap_pts = test_pts(id_overlap,:);
        % create the grid of overlapping points to use in the gridded
        % interpolant
        x = unique(overlap_pts(:,1));
        y = unique(overlap_pts(:,2));
        [X, Y] = ndgrid(x,y);
        % evaluate the gridded interpolant of the lower box using the
        % points that are within the test box
        lower_vals = fh_lower.F(X,Y);
        % find the test points within the test box and compare with the
        % newly found points
        [~,idx,~] = intersect(x_test,x);
        [~,idy,~] = intersect(y_test,y);
        % int_vals is a dummy variable to make the indexing of the boxes
        % the same for comparison
        int_vals = test_vals(idx,idy);
        % replace all test values that are greater than the lower values
        int_vals(int_vals>lower_vals) = lower_vals(int_vals>lower_vals);
        % place the int_vals back into the test values
        test_vals(idx,idy) = int_vals;
    end
    % put the new points back in the edgefunction
    fh{kk}.F = griddedInterpolant(X_test,Y_test,test_vals);
    
    
    
end
