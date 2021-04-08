function [g10,g11,h11]=mc_igrf(dateNums)
% MC_IGRF Returns value of dipole coefficients for IGRF at given dates and
% times.
% 
% IGRF values are read from file 'igrf.mat' if present, if not this file
% will be created from the latest file of type 'igrf<ver.>coeffs.txt'
% which is currently available from
%    https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
% 
% dateNums  Vector of dates and times to load coefficients for. Either
%           MATLAB serial date numbers or cell of strings that can be
%           understood by DATEVEC e.g. {'2015-01-01 12:42:30', '2016-05-12
%           00:45:11'}
% g10       Vector of g10 coefficients at dateNums, in nanoTesla
% g11       Vector of g11 coefficients at dateNums, in nanoTesla
% h11       Vector of h11 coefficients at dateNums, in nanoTesla

% When a new IGRF coefficient file is available, place it in this
% directory, delete 'igrf.mat', and a new 'igrf.mat' will be generated on
% first call to this function.

% References:
%
% Alken et al, International Geomagnetic Reference Field: the 13th
% generation, Earth Planets Space, 2020.
% 
% https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html

% Will Brown, January 2020
%
% Jan/2020 - brought this function directly into m_map without changes.


if ~exist('igrf.mat', 'file') % Create .mat file if not present
    currLoc=mfilename('fullpath');
    % Use most recent file version number if multiple present
    fileList=dir([currLoc(1:end-7),'igrf*coeffs.txt']);
    [refTimes,igrfg10,igrfg11,igrfh11]=readIGRFFile( ...
        fullfile(fileList(end).folder,fileList(end).name));
    save([currLoc(1:end-7),'igrf.mat'],'refTimes','igrfg10','igrfg11', ...
        'igrfh11')
else
    load('igrf.mat') % vars: refTimes, igrfg10, igrfg11, igrfh11
end

% Convert dateNums to Matlab serial DATENUM if char strings given
if ischar(dateNums) || iscell(dateNums)
    dateNums = datenum(dateNums);
end
% Convert Matlab serial DATENUM to decimal years
nDateTimes = size(dateNums,1);
dateVecs = datevec(dateNums);
decYears = dateVecs(:,1)+(dateNums-datenum(...
    [dateVecs(:,1),ones(nDateTimes,1),ones(nDateTimes,1)]))...
    ./(365+double((~mod(dateVecs(:,1),4) & mod(dateVecs(1,1),100)) ...
    | (~mod(dateVecs(:,1),400)))); % with leap years

% Check validity of model vs input times
if min(decYears) < refTimes(1) || max(decYears) > refTimes(end)
    error(['mc_igrf: IGRF version loaded is only valid between ', ...
        num2str(refTimes(1)),'-01-01 00:00:00 and ', ...
        num2str(refTimes(end)),'-01-01 00:00:00.']);
end

g10 = interp1(refTimes,igrfg10,decYears);
g11 = interp1(refTimes,igrfg11,decYears);
h11 = interp1(refTimes,igrfh11,decYears);

end % funnction MC_IGRF

function [refTimes,igrfg10,igrfg11,igrfh11]=readIGRFFile(igrfFile)
% READIGRFFILE Read IGRF coefficient file and returns times and values of
%  dipole coefficients, including end of 5 year predictive component.
% 
% igrfFile  IGRF coefficient file of type 'igrf<ver.>coeffs.txt'
% refTimes  Vector of coefficient times, decimal years
% igrfg10   Vector of g10 coefficients at refTimes, nanoTesla
% igrfg11   Vector of g11 coefficients at refTimes, nanoTesla
% igrfh11   Vector of h11 coefficients at refTimes, nanoTesla

% Assuming official 'igrf<ver.>coeffs.txt' file format goes unchanged in
% future...
igrf=importdata(igrfFile, ' ', 4);

% Pull dates list from header text
refTimes=igrf.textdata(4,4:end-1);
refTimes=sprintf('%s ',refTimes{:});
refTimes=sscanf(refTimes,'%f');
refTimes=[refTimes;refTimes(end)+5]; % Include 5year predictive period

% Get main field dipole coeffs
igrfg10=igrf.data(1,3:end-1)';
igrfg11=igrf.data(2,3:end-1)';
igrfh11=igrf.data(3,3:end-1)';

% Add secular variation dipole coeffs for final predictive 5year period
igrfg10=[igrfg10;igrfg10(end)+5*igrf.data(1,end)];
igrfg11=[igrfg11;igrfg11(end)+5*igrf.data(2,end)];
igrfh11=[igrfh11;igrfh11(end)+5*igrf.data(3,end)];

% Check equal numbers of times and coeffs found
if ~(numel(refTimes) == numel(igrfg10) ...
        && numel(igrfg10) == numel(igrfg11) ...
        && numel(igrfg11) == numel(igrfh11))
    error(['mc_igrf: readIGRFFile: Unequal numbers of times and'; ...
        ' coefficients read from %s.'],igrfFile)
end

end % function READIGRFFILE