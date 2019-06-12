function merge = mergeTiles(tiles,globaltile,inputname,checktimestep_params) 

global MAP_PROJECTION MAP_COORDS MAP_VAR_LIST

if nargin == 1
    % Merge mesh tiles together 
    % complex looping strategy is faster than the trivial addition because
    % it has fewer additions with a large mesh
    % note that plus can now handle connected and non-connected tiles
    
    numTiles = length(tiles);
    n = 0;
    for i = 1: 2 : numTiles
        dmy = load(tiles{i}); 
        dmy = dmy.m1;
        n = n + 1;
        disp(['merge #' num2str(n)])
        merge{n} = dmy;
        if i + 1 <= numTiles
            dmy = load(tiles{i+1}) ; 
            dmy = dmy.m1; 
            merge{n} = dmy + merge{n};
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
    
else
    if nargin < 3 || isempty(inputname)
       inputname = {'m' 'm'}; 
    end
    if nargin == 4
        DT = checktimestep_params(1); CFL = checktimestep_params(2);
    end
    % Merge mesh tiles together with a global one
    load(globaltile); 
    merge = eval(inputname{1});
    if exist('ii'); is = ii+1; else is = 1; end
    for ii = is:length(tiles)
        load(tiles{ii})
        m = eval(inputname{2}); buf = 0.5; clean_cutoff = 1e-4;
        bbox = [min(m.p)-buf; max(m.p)+buf]';
        % extract the regional domain from global domain (and get the inverse)
        mregion = ExtractSubDomain(merge,bbox);
        setProj(mregion,1,merge.proj.name) ;
        mregion.proj    = MAP_PROJECTION ;
        mregion.coord   = MAP_COORDS ;
        mregion.mapvar  = MAP_VAR_LIST ;   
        mregion.bd = []; mregion.op = [];
        mregion.f13 = []; mregion.f15 = []; mregion.f24 = [];
        mregionold = mregion; 
        mregion = clean(mregion,0,0,0,clean_cutoff,[],[],[],0);
        if size(mregionold.p,1) ~= size(mregion.p,1)
            nearest = ourKNNsearch(mregionold.p',mregion.p',1);
            mregion.b = mregionold.b(nearest);
            if ~isempty(mregion.bx)
                mregion.bx = mregionold.bx(nearest);
                mregion.by = mregionold.by(nearest);
            end
        end
        % add together the regional domain to the high resolution domain
        mcom = plus(m,mregion,0); 
        % make sure timestep condition is met for the combined mesh
        if nargin == 3
            mcom = CheckTimestep(mcom,DT,CFL);
        end
        % trivially transfer slope to combined domain 
        nearest = ourKNNsearch(mregion.p',mcom.p',1);
        if ~isempty(mregion.bx)
            mcom.bx = mregion.bx(nearest);
            mcom.by = mregion.by(nearest);
        end
        mglobal = ExtractSubDomain(merge,bbox,1);
        % trivially merge the combined region and cutout global domain
        merge = xor(mcom,mglobal);
        
        if ~mod(ii,10) 
           save('temp_merge','merge','ii','tiles');
        end
    end    
    
end

end
