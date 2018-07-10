function [idx, dst] = ourKNNsearch(dataset,testset,k) 
% This wrapper is used because we cannot pass an ANN class object to
% parfeval in MATLAB R2017A. 
% See reference (https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ANNmanual_1.1.pdf) returns the squared distances 
% ANN: A Library for
% Approximate Nearest Neighbor Searching
% David M. Mount and Sunil Arya
% Version 1.1.2
% Release Date: Jan 27, 2010
anno = ann(dataset); 
[idx,~] = ksearch(anno, testset,k,0); 
idx = idx'; 
datax = dataset(1,:)';
datay = dataset(2,:)';
testx = repmat(testset(1,:)',1,k);
testy = repmat(testset(2,:)',1,k);
dst = sqrt((testx-datax(idx)).^2 + (testy-datay(idx)).^2);
anno = close(anno); 
end