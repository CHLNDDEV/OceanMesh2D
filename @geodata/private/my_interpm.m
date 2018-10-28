function [latout,lonout] = my_interpm(lat,lon,maxdiff)
% This function  fills in any gaps in latitude (lat) or longitude (lon) data vectors 
% that are greater than a defined tolerance maxdiff apart in either dimension.
% lat and lon should be row vectors.
%
% latout and lonout are the new latitude and longitude data vectors, in which any gaps
% larger than maxdiff in the original vectors have been filled with additional points.
%
% by Jindong Wang on March 1, 2018
%
if isrow(lat)
    lat=lat';
end
if isrow(lon)
    lon=lon';
end
ny=size(lat,1);
nx=size(lon,1);
if nx~=ny || nx<2 || ny<2 || maxdiff<10^-8
    error('Error: Wrong input!');
end

dlat=abs(lat(2:end)-lat(1:end-1));
dlon=abs(lon(2:end)-lon(1:end-1));
nin=ceil(max([dlat,dlon],[],2)/maxdiff)-1;
sumnin=sum(nin,'omitnan');
if sumnin==0
    disp('No insertion needed.');
    latout=lat;
    lonout=lon;
    return
end
nout=sumnin+nx;
latout=nan(nout,1);
lonout=nan(nout,1);

n=1;
for i=1:nx-1;
    ni=nin(i);
    if ni==0 || isnan(ni)
        latout(n)=lat(i);
        lonout(n)=lon(i);
        nstep=1;
    else
        ilat=linspace(lat(i),lat(i+1),ni+2);
        ilon=linspace(lon(i),lon(i+1),ni+2);
        latout(n:n+ni)=ilat(1:ni+1);
        lonout(n:n+ni)=ilon(1:ni+1);
        nstep=ni+1;
    end
    n=n+nstep;
end
latout(end)=lat(end);
lonout(end)=lon(end);
        