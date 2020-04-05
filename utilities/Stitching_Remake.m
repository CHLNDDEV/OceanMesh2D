function poly = Stitching_Remake(m,upc,dir,islands)
    % Stiches the shoreline from an ocean only msh object (m) to an upper
    % contour (gdat) inside the bbox to get a Nan-delimited polygon for
    % meshing the floodplain
    
    % dir is a direction may need to try and test (0 or 1)
    % we assume that everything is going clockwise

    % get the ocean only shoreline
    cell2 = extdom_polygon(extdom_edges2(m.t,m.p),m.p,-1,0);
    ocean_shoreline = cell2mat(cell2');
    
    % save the first part of ocean shoreline 
    % which we stitch with to the upper contour
    Inan = find(isnan(ocean_shoreline(:,1)),1,'first');
    ocean_first = ocean_shoreline(1:Inan-1,:);
    % save the rest (islands)
    ocean_shoreline(1:Inan-1,:) = [];
    
    % get nearest point of ocean_first to upc(end)
    d = sum((ocean_first - upc(end,:)).^2,2);
    [~,Ie] = min(d);
    % removing repeated first and last point in process of rearranging
    ocean_first = [ocean_first(Ie:end-1,:); ocean_first(1:Ie-1,:)];
    
    % get nearest point of ocean_first to upc(1)
    d = sum((ocean_first - upc(1,:)).^2,2);
    [~,If] = min(d);
    
    % construct the polygon
    if dir 
        poly = [upc; ocean_first(1:If,:); upc(1,:)]; 
    else
        poly = [upc; ocean_first(1,:); ocean_first(end:-1:If,:); upc(1,:)]; 
    end
    
    if islands
        % return the full polygon
        poly = [poly; NaN NaN; ocean_shoreline];
    end
end