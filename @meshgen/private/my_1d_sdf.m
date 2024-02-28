function [dist] = my_1d_sdf(p,A,B)
distA = abs(p - A);
distB = abs(p - B);
dist = -min(distA,distB);
if p < A
    dist = -1*dist;
elseif p > B
    dist = -1*dist;
end
end