f14 = 'fort.14';
f24 = 'fort.24';
lonlat0 = [  0 0];
saldata = 'FES2014' ;
avisoloc = '..';
f15tipname = {'Q1','O1','P1','K1','N2','M2','S2','K2'};
tic
FnGlobal_SAL_to_fort24( f24, f14, f15tipname, lonlat0, avisoloc, saldata )
toc