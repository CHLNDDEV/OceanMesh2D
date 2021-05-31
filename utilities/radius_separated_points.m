function ind = radius_separated_points(long,lat,radius)
%  ind = radius_separated_points(long,lat,radius)
%  get the indices of points in [long, lat] to keep 
%  so that they are at least separated by a distance of
%  user-defined 'radius'

% computing rangesearch (ourRangeSearch significantly slower due to for loop)
if ~exist('rangesearch')
    K = ourRangeSearch([long lat]',[long lat]',radius);
else
    K = rangesearch([long lat],[long lat],radius);
end
KL = cellfun(@length,K);
rr = find(KL > 1);
% loop until all points have closest point that is only itself
del = [];
while ~isempty(rr)
   rI = randi(length(rr),1);
   kk = K{rr(rI)};
   % save the locations to delete
   del = [del kk(2:end)];
   KL(kk) = 1; % change the KL of kk to 1
   rr = find(KL > 1);
end
del = unique(del);
ind = setdiff(1:length(long),del);
 
