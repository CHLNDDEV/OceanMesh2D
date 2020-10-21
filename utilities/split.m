% Function for adding points to an untriangulated set of points based on
% nearest distance versus the desired edgefx 
function [new_points]=split(points,fh)
    % find the point closest to each point (that isn't the *same* point)
    [idx, ~] = ourKNNsearch(points',points',2);
    % the ideal spacing between points
    ideal_dst = fh(points);
    % where the dst is 2*ideal_dist, add a point in between
    long   = zeros(length(points)*2,1);
    lat    = zeros(length(points)*2,1);
    long(1:2:end) = points(idx(:,1),1);
    long(2:2:end) = points(idx(:,2),1);
    lat(1:2:end)  = points(idx(:,1),2);
    lat(2:2:end)  = points(idx(:,2),2);
    L = m_lldist(long,lat);
    L = L(1:2:end)*1e3;            % L = lengths in meters
    splits = L > 1.5*ideal_dst;
    disp(['Number of splits near multiscale nest: ' num2str(sum(splits))])
    new_points = (points(idx(splits,1),:) + points(idx(splits,2),:))./2.0;
end


