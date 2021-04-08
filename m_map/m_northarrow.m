function h=m_northarrow(long,lat,scale,varargin)
% M_NORTHARROW draws a north arrow on a map
%   M_NORTHARROW(LON,LAT,SCALE) draws a northward-pointing
%   arrow on a map. THe arrow is located at the single geographic
%   coordinates represented by scalar (LON,LAT), and extends across 
%   SCALE degrees of latitude.
%
%   M_NORTHARROW(LON,LAT,SCALE,...PARAM,VALUE...) lets you
%   specify patch parameters (e.g., 'facecolor', or 'edgecolor')
%
%   Some special parameters affect the look of the arrow:
%       'type' : range 1-4, chooses different designs (default 1)
%       'aspect': If >1, makes objects thinner (above 3 no recommended)
%                 If <1, makes object fatter (below .6 not recommended)
%

% R. Pawlowicz Jan/2020


global MAP_PROJECTION  

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

% Parse inputs for arguments

% Defaults
type=4;
aspect=1;

k=1;
while k<length(varargin)
    switch varargin{k}(1:3)
       case 'typ'
          type=varargin{k+1};
          varargin([k k+1])=[];
        case 'asp'
          aspect=varargin{k+1};
          varargin([k k+1])=[];
        otherwise  % Pass arguments through 
          k=k+2;
    end
end


% Get the northward direction
[x,y]=m_ll2xy([long long],[lat lat+.01],'clip','off');
% base vector
dN=(diff(x)+1i*diff(y))*100;

 
gry=[.7 .7 .7];

switch type
    case 1   % Black no-stem arrow

        % The arrows
        za=1i*[0.1 0.1 1.5 0 -1.5 -.1 -.1 .1]/aspect+[-.5 1 0 5 0  1 -.5 -.5]  - 2.5;
        % The N
        zN=1i*[-.7 -.7 .7 .7 .7 -.7 -.7]+[8 6 8 6 8 6 8] - 3;

        % Normalize it
        zN =(zN )/5.*dN*scale;
        zab=(za )/5.*dN*scale;
        zaw=NaN;
        
    case 2   % Black and white no-stem arrow
                % The arrows
        za=1i*[  0 1.5 0 0]/aspect+[ 1 0 5  1]  - 2.5;
         % The N
        zN=1i*[-.7 -.7 .7 .7 .7 -.7 -.7]+[8 6 8 6 8 6 8] - 3;

        % Normalize it
        zN =(zN )/5.*dN*scale;
        zaw=(za )/5.*dN*scale;
        zab=(za')/5.*dN*scale;

    case 3    % 4-way rose
        
        za=[0 1/aspect 6 0]+1i*[0 1 0 0]/aspect;
        za=[za 1i*za -za -1i*za];
          % The N
        zN=1i*[-.7 -.7 .7 .7 .7 -.7 -.7]+[8 6 8 6 8 6 8] ;
        % Normalize it
        zN =(zN )/12.*dN*scale;
        zaw=(za )/12.*dN*scale;
        zab=(za')/12.*dN*scale;

    case 4     % 8-way rose
        sfac=2/aspect;
        za=[0 sfac*cosd(22.5) 6 0 sfac*cosd(67.5) 3 0]+1i*[0 sfac*sind(22.5) 0 0 sfac*sind(67.5) 3 0];
        za=[za 1i*za -za -1i*za];
          % The N
        zN=1i*[-.7 -.7 .7 .7 .7 -.7 -.7]+[8 6 8 6 8 6 8] ;
        % Normalize it
        zN =(zN )/12.*dN*scale;
        zaw=(za )/12.*dN*scale;
        zab=(za')/12.*dN*scale;
             
end

h(1)=patch(x(1)+real(zN), y(1)+imag(zN),'k', 'linewi',2,'edgecolor','k',varargin{:},'clip','off');
h(2)=patch(x(1)+real(zab),y(1)+imag(zab),gry,'linewi',1,'edgecolor','k',varargin{:},'clip','off');
h(3)=patch(x(1)+real(zaw),y(1)+imag(zaw),gry,'linewi',1,'edgecolor','k',varargin{:},'facecolor','w','clip','off');
  
 
