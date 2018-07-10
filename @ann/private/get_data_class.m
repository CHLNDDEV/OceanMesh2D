function [ccls cfun] = get_data_class()
% get class supported by ann 
d = annmex();
ccls = class(d);

cfun = @(x) cast(x, ccls);
