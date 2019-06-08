function merge = mergeTiles(tiles) 
    % Merge mesh tiles together 
    % complex looping strategy is faster than the trivial addition because
    % it has fewer additions with a large mesh
    % note that plus can now handle connected and non-connected tiles
    
    numTiles = length(tiles);
    n = 0;
    for i = 1: 2 : numTiles
        dmy = load(tiles{i});
        n = n + 1;
        disp(['merge #' num2str(n)])
        merge{n} = dmy.m;
        if i + 1 <= numTiles
            dmy = load(tiles{i+1}) ;
            merge{n} = dmy.m + merge{n};
        end
    end
    
    while length(merge) > 1
        n = 0;
        for i = 1 : 2 : length(merge)
            n = n + 1;
            disp(['merge #' num2str(n)])
            if i + 1 <= length(merge) 
                merge{n} = merge{i+1} + merge{i};
            else
                merge{n} = merge{i}; 
            end
        end 
        merge(n+1:end) = [];
    end
    merge = merge{1};
end
