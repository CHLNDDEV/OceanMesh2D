function [idx, dst] = WrapperForKsearch(dataset,testset) 
% This wrapper is used because we cannot pass an ANN class object to
% parfeval in MATLAB R2017A. 
% See reference (https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ANNmanual_1.1.pdf) returns the squared distances 
% ANN: A Library for
% Approximate Nearest Neighbor Searching
% David M. Mount and Sunil Arya
% Version 1.1.2
% Release Date: Jan 27, 2010
anno = ann(dataset); 
[idx,~] = ksearch(anno, testset,1,0); 
idx = idx'; 
dst = sqrt((testset(1,:)'-dataset(1,idx)').^2 + (testset(2,:)'-dataset(2,idx)').^2);
anno = close(anno); 
end