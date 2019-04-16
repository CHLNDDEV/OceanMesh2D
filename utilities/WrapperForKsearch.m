function [idx, dst] = WrapperForKsearch(anno,dataset,testset,project) 
% This wrapper is used because we cannot pass an ANN class object to
% parfeval in MATLAB R2017A. 
% See reference (https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ANNmanual_1.1.pdf) returns the squared distances 
% ANN: A Library for
% Approximate Nearest Neighbor Searching
% David M. Mount and Sunil Arya
% Version 1.1.2
% Release Date: Jan 27, 2010

% Do some projection
if nargin < 4
   project = 0; 
end
if project
%      [dataset(:,1),dataset(:,2)] = m_ll2xy(dataset(:,1),dataset(:,2));
       [testset(:,1),testset(:,2)] = m_ll2xy(testset(:,1),testset(:,2));
%      dataset(isnan(dataset(:,1)),:) = [];
% %     % This line removes the line that can appear in the center for
% %     % steoreographic projection from the bbox
%      dataset(abs(dataset(:,1)) < 1e-6,:) = [];
 end
%anno = ann(dataset'); 
[idx,~] = ksearch(anno, testset',1,0); 
idx = idx'; 
dst = sqrt((testset(:,1)-dataset(idx,1)).^2 + (testset(:,2)-dataset(idx,2)).^2);
%anno = close(anno); 
end