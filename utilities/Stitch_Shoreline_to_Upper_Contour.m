function poly = Stitch_Shoreline_to_Upper_Contour(m,gdat,bbox)
    % Stiches the shoreline from an ocean only msh object (m) to an upper
    % contour (gdat) inside the bbox to get a Nan-delimited polygon for
    % meshing the floodplain

    if nargin < 3
        bbox = [-inf inf; -inf inf];
    end

    % get the floodplain inside bbox
    outer_contour = gdat.outer;
%     % remove the bounding box
%     Inan = find(isnan(outer_contour(:,1)),1,'first');
%     outer_contour(Inan:end,:) = [];
%     outer_contour(1:Inan,:) = [];
%     del = outer_contour(:,1) < bbox(1,1) | outer_contour(:,1) > bbox(1,2) | ...
%           outer_contour(:,2) < bbox(2,1) | outer_contour(:,2) > bbox(2,2);
%     outer_contour(del,:) = [];   
%     outer_contour(isnan(outer_contour(:,1)),:) = [];
%     % reorder if necessary
%     [mo1,I] = max(sum(diff(outer_contour).^2,2));
% 	outer_contour_test = [outer_contour(I+1:end,:); outer_contour(1:I,:)];
%     [mo2,~] = max(sum(diff(outer_contour_test).^2,2));
%     if mo2 < mo1
%         outer_contour = outer_contour_test;
%     end
    %K = convhull(outer_contour(:,1),outer_contour(:,2));
    %polycheck = outer_contour(K,:);
    polycheck = outer_contour;
    
    % get the ocean only shoreline
    cell2 = extdom_polygon(extdom_edges2(m.t,m.p),m.p,-1,0);
    ocean_shoreline = cell2mat(cell2');
    
    % save the first part of ocean shoreline 
    % which we stitch with to the upper contour
    Inan = find(isnan(ocean_shoreline(:,1)),1,'first');
    ocean_first = ocean_shoreline(1:Inan-1,:);
    
    % save the rest (islands)
    ocean_shoreline(1:Inan-1,:) = [];
    % remove unneccessary islands that are outside of bbox
    Inans = find(isnan(ocean_shoreline(:,1))); del = 1;
    for ii = 1:length(Inans)-1
       %if all(ocean_shoreline(Inans(ii)+1:Inans(ii+1)-1,1) < bbox(1,1) | ...
       %       ocean_shoreline(Inans(ii)+1:Inans(ii+1)-1,1) > bbox(1,2) | ...
       %       ocean_shoreline(Inans(ii)+1:Inans(ii+1)-1,2) < bbox(2,1) | ...
       %       ocean_shoreline(Inans(ii)+1:Inans(ii+1)-1,2) > bbox(2,2))
       if all(~inpoly(ocean_shoreline(Inans(ii)+1:Inans(ii+1)-1,:),polycheck))
          del = [del Inans(ii)+1:Inans(ii+1)];
       end
    end
    ocean_shoreline(del,:) = [];
    
    % get the main shoreline inside bbox
    %del = ocean_first(:,1) < bbox(1,1) | ocean_first(:,1) > bbox(1,2) | ...
    %      ocean_first(:,2) < bbox(2,1) | ocean_first(:,2) > bbox(2,2);
    del = ~inpoly(ocean_first,polycheck);
    ocean_first(del,:) = [];
    % reorder if necessary
    [mo1,I] = max(sum(diff(ocean_first).^2,2));
	ocean_first_test = [ocean_first(I+1:end,:); ocean_first(1:I,:)];
    [mo2,~] = max(sum(diff(ocean_first_test).^2,2));
    if mo2 < mo1
        ocean_first = ocean_first_test;
    end
    
    % stitch ocean_first to outer_contour 
    poly = ocean_first;
    OFEOCF = ocean_first(end,:) - outer_contour(1,:);
    OFEOCE = ocean_first(end,:) - outer_contour(end,:);
    if hypot(OFEOCF(1),OFEOCF(2)) <= hypot(OFEOCE(1),OFEOCE(2))
        poly = [poly; outer_contour];
    else
        poly = [poly; flipud(outer_contour)];
    end
    % remove any unique points
    poly = unique(poly,'rows','stable');
    
    % make sure first point equals last
    poly(end+1,:) = poly(1,:);
    
    % return the full polygon
    poly = [poly; NaN NaN; ocean_shoreline];
end