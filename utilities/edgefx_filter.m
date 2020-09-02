function fh = edgefx_filter(fh,gdat)
% This function sorts through layered edgefunctions and puts the minimum
% target edgelength in the topmost edgefunction. This ensures that when
% the mesh is generated, the minimum target edgelength from all boxes is
% used.
% Coleman Blakely
% January 22, 2020
if length(fh) == 1
    disp('No need to do this silly')
end
% preallocate polyshape vector for finding whether boundaries are
% overlapping
S = cell(length(fh)-1,1);
% build polyshapes starting with last (top) edgefunction and working down
for kk = length(fh):-1:2
    S{kk} = polyshape(fh{kk}.boubox);
end
% Loop through edgefunctions to 1) check if they overlap and 2) if they
% overlap find the points that lie within the upper edgefunction and
% compare them to the ones below
for kk = length(fh):-1:2
    fh_test = fh{kk};
    % make the points to use with inpoly to find the overlapping points
    [test_pts, ~] = edgefx2pts(fh_test);
    test_vals = fh_test.F.Values;
    x_test = fh_test.F.GridVectors{1};
    y_test = fh_test.F.GridVectors{2};
    [X_test, Y_test] = ndgrid(x_test,y_test);
    % Loop through all boxes that may be below the test box
    for ii = kk-1:-1:2
        % check if the test box overlaps the below box
        if ~overlaps(S{kk},S{ii})
            continue
        end
        fh_lower = fh{ii};
        % find the test points that are within the lower boundary
        id_overlap = inpoly(test_pts,fh_lower.boubox(1:end-1,:));
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
% now resmooth the edgefunctions using the defined grade to eliminate
% any discontinuities, etc.
% most of this was pulled from stuff that William and Keith did. I do
% not know exactly how all of it works but it appears to do the trick.
% CPB
for kk = length(fh):-1:1
    obj = fh{kk};
    feat = gdat{kk};
    hh_m = obj.F.Values;
    [xg, yg] = CreateStructGrid(obj);
    if all(abs(obj.bbox(1,:)) == 180)
        % for global mesh make it cyclical
        hh_m = [hh_m(end,:); hh_m; hh_m(1,:)];
    end
    hfun = reshape(hh_m',[numel(hh_m),1]);
    
    dx = obj.h0*cosd(min(yg(1,:),85)); % for gradient function
    dy = obj.h0;                       % for gradient function
    
    % make g a function of space
    dmy     = xg*0 ;
    for param = obj.g'
        if numel(param)==1 && param~=0
            lim   = obj.g(1);
            dmy  = dmy + lim ;
        else
            lim  = param(1);
            dp1  = param(2);
            dp2  = param(3);
            
            limidx = (feat.Fb(xg,yg) < dp1 & ...
                feat.Fb(xg,yg) > dp2) ;
            
            dmy( limidx ) = lim;
        end
    end
    if all(abs(obj.bbox(1,:)) == 180)
        % for global mesh make it cyclical
        dmy = [dmy(end,:); dmy; dmy(1,:)];
    end
    fdfdx = reshape(dmy',[numel(dmy),1]);
    clearvars dmy;
    [hfun,flag] = limgradStruct(obj.ny,dx,dy,hfun,...
        fdfdx,sqrt(length(hfun)));
    if flag == 1
        disp('Gradient relaxing converged!');
    else
        error(['FATAL: Gradient relaxing did not converge, '
            'please check your edge functions']);
    end
    % reshape it back
    hh_m = reshape(hfun,obj.ny,[])';
    clearvars hfun fdfdx
    if all(abs(obj.bbox(1,:)) == 180)
        % for global mesh make it cyclical
        hh_m = hh_m(2:end-1,:);
        hh_m(end,:) = hh_m(1,:);
    end
    obj.F = griddedInterpolant(xg,yg,hh_m,'linear','nearest');
    
    clearvars xg yg
    fh{kk} = obj;
end
end

%% This is a quick little function to turn an edgefuncton into points and values that can be used with fastscatter
function [pts, vals] = edgefx2pts(fh)
% takes an edgefunction from OceanMesh2D and outputs a series of points and
% minimum element sizes for each point
x = fh.F.GridVectors{1};
y = fh.F.GridVectors{2}';
val = fh.F.Values;
vals = reshape(val,[],1);
[lat, lon] = meshgrid(x,y);
c = cat(2,lat',lon');
pts = reshape(c,[],2);
end