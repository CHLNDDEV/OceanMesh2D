function idx = ourRangeSearch(dataset,testset,r,max_nei) 
% idx = ourRangeSearch(dataset,testset,r,max_nei) 
% dataset/testset = d (dimension) x N (length) arrays
% r - fixed radius to search
% max_nei - optional setting for maximum allowable neighbors in radius for
% each query point. Default will be 5% of lenght of dataset if not set

if size(testset,1) ~= size(dataset,1)
    error('Dimensions of dataset and testset (first axis) are not equal')
end
query_pts = size(testset,2);
if nargin < 4
    max_nei = floor(size(dataset,2)*0.05); % try 5% of length ofdataset
end

anno = ann(dataset); 
idx = cell(query_pts,1);
% loop over all points and add to cell
for pts = 1:query_pts
    idx{pts} = frsearch(anno, testset(:,pts), r, max_nei,0,1);
end
if any(cellfun(@length,idx) == max_nei)
    warning('Number of neighbours in radius exceeded limit for at least one query point')
end
anno = close(anno); 
end

