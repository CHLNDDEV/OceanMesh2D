function [X,Y,I]=m_ll2xy(varargin)
% M_LL2XY Converts long,lat to X,Y coordinates using the current projection
%         [X,Y]=m_ll2xy(LONGITUDE,LATITUDE);
%
%         Extra properties can be added after the latitude variable:
%          ...,'clip', ( 'on' | 'off' | 'patch'  | 'point' )
%         where normal clipping sets out-of-range values to NaN, and patch
%         clipping sets out-of-range values to border values (useful for
%         drawing patches). A border point is interpolated for line
%         segments crossing the border; line segments are assumed for
%         one-dimensional input arrays and for the columns of 2-dimensional
%         arrays. The interpolation can be disabled using the 'point'
%         option (i.e. data is composed of discrete unrelated points).
%
%         [X,Y,I]=m_ll2xy(...) returns an index to tell you which
%         points have been modified (I=1).

% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 4/DEc/11 - isstr to ischar
% 1/Nov/12 - added another geomagnetic coordinate system
% 10/Jan/20 - handling for any geomagnetic coord system

global MAP_PROJECTION MAP_COORDS


if isempty(MAP_PROJECTION)
   disp('No Map Projection initialized - call M_PROJ first!');
   return;
end


if nargin==0 || ischar(varargin{1})
  disp(' Usage');
  disp(' [X,Y]=m_ll2xy(LONGITUDES,LATITUDES <,''clip'',( ''on''|''off''|''patch'' | ''point'' ) >)');
else
  if strcmp(MAP_COORDS.name.name,MAP_PROJECTION.coordsystem.name.name) && ...  % using same coord system as projection
          MAP_COORDS.name.mdate==MAP_PROJECTION.coordsystem.name.mdate
     % Sneaky way of making default clipping on (sneaky 'cause only the 4th
     % input parameter is checked for the clipping property)
     [X,Y,I]=feval(MAP_PROJECTION.routine,'ll2xy',varargin{:},'clip','on');
     
  elseif strcmp(MAP_COORDS.name.name,'geographic')    % Projection uses geomagnetic coords, data uses geographic
     [LONG,LAT]=mc_coords('XYgeo2mag',varargin{1:2});
  %   plot(varargin{1},varargin{2},'.-',LONG,LAT,'.-');
     % The problem with this conversion is that the converted coordinates
     % sometimes jump by +/- 360 degrees. Here I try to reduce that problem
     % by making lines 'continuous' without those jumps.
     % Right now I am only dealing with VECTOR inputs.
     Lin=varargin{1};
     if isvector(Lin)
         kk=find(isnan(Lin));
         if any(kk)   % Lines separated by NaNs
             for k=1:length(kk)-1
                k2=kk(k)+1:kk(k+1)-1;
                if any(k2)
                    if Lin(k2(1))>180+360
                         ioff=360*2;
                    elseif Lin(k2(1))>180
                         ioff=360;
                     elseif Lin(k2(1))<-180
                         ioff=-360;
                     elseif Lin(k2(1))<-180-360
                         ioff=-360*2;
                     else
                         ioff=0;
                    end
                  %  fprintf('%d ',ioff);
                     if iscolumn(Lin)
                         LONG(k2)=LONG(k2)+ioff+...
                             ([0;cumsum(diff(LONG(k2))<-250)-cumsum(diff(LONG(k2))>250)])*360;
                     else
                         LONG(k2)=LONG(k2)+ioff+...
                             ([0,cumsum(diff(LONG(k2))<-250)-cumsum(diff(LONG(k2))>250)])*360;
                     end
                end
             end
         else   % No NaNs separating lines
           %  Lin(1)
             if Lin(1)>180+360
                 ioff=360*2;
             elseif Lin(1)>180
                 ioff=360;
             elseif Lin(1)<-180
                 ioff=-360;
             elseif Lin(1)<-180-360
                 ioff=-360*2;
             else
                 ioff=0;
             end
          %   fprintf('%d ',ioff);
             if iscolumn(Lin)
                 LONG=LONG+ioff+...
                     ([0;cumsum(diff(LONG)<-250)-cumsum(diff(LONG)>250)])*360;
             else
                 LONG=LONG+ioff+...
                     ([0,cumsum(diff(LONG)<-250)-cumsum(diff(LONG)>250)])*360;
             end
         end
     end  
 %    line(LONG,LAT,'color','r','linewi',2);
     args={varargin{3:end},'clip','on'};
     [X,Y,I]=feval(MAP_PROJECTION.routine,'ll2xy',LONG,LAT,args{:});
     
  else      % Data uses geomagnetic coords, map uses geographic 
     [LONG,LAT]=mc_coords('XYmag2geo',varargin{1:2});
     args={varargin{3:end},'clip','on'};
     [X,Y,I]=feval(MAP_PROJECTION.routine,'ll2xy',LONG,LAT,args{:});
     
  end 
end

if nargout==0
 clear X Y I
end

