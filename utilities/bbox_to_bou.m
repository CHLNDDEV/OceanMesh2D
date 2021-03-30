function bou = bbox_to_bou(bbox) 
   % bou = bbox_to_bou(bbox) 
   % Converts a bbox = [lon_min lon_max;
   %                    lat_min lat_max];
   % to a polygon with NaNs at the end

   bou = [bbox(1,1) bbox(2,1);
          bbox(1,1) bbox(2,2); 
          bbox(1,2) bbox(2,2);
          bbox(1,2) bbox(2,1); 
          bbox(1,1) bbox(2,1); 
          NaN       NaN];
end
