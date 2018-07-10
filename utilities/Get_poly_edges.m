function edges = Get_poly_edges( poly )
%   edges = Get_poly_edges( poly )
%   Finds the edges for input to inpoly which is required for multiple 
%   polygons stitched together using nans
    shpEnd = find(isnan(poly(:,1))); 
    shpEnd = vertcat(0,shpEnd); 
    edges = nan(length(poly(:,1))-length(shpEnd),2); 
    count = 1; 
    for j=1:length(shpEnd)-1 
        endCount = count+length((shpEnd(j)+1:shpEnd(j+1)-2)); 
        edges(count:endCount,:) = [(shpEnd(j)+1:shpEnd(j+1)-2)' ... 
        (shpEnd(j)+2:shpEnd(j+1)-1)';shpEnd(j+1)-1 shpEnd(j)+1]; 
        count = endCount+1; 
    end

end

