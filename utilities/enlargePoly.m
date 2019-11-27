function [xout,yout]=enlargePoly(x,y,frac)
% Enlarges a NaN-delimited polygon by a fractional amount
xout = frac*x+(1-frac)*nanmean(x);
yout = frac*y+(1-frac)*nanmean(y);

