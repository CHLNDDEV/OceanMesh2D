function [longOUT,latOUT,phiVecOUT,thetaVecOUT]=mc_coords(optn,longIN,latIN,phiVecIN,thetaVecIN)
% MP_COORD Converts between coordinate systems based on different
%  poles. Generally used in space physics to plot things in
%  geomagnetic (dipole) coordinate systems.
%
%  This function should not be used directly; instead it is
%  is accessed by various high-level functions named M_*.
%
%  Vector rotations can be carried out using the following:
%
%   [latOUT,longOUT,phiVecOUT,thetaVecOUT]=M_GEO2MAG(latIN,longIN,phiVecIN,thetaVecIN)
%
%  where
%
% thetaVecIN - north component of the vector in geographic coordinates
% phiVecIN   - east component of the vector in geographic coordinates
% thetaVecOUT - north component of the vector in geomagnetic coordinates
% phiVecOUT   - east component of the vector in geomagnetic coordinates
%

% References:
%
% Hapgood, M.A., Space Physics Coordinate Transformations:
% A User Guide, Planet. Space Sci., Vol. 40, N0. 5, 1992.

% Antti Pulkkinen, September 2001.
% R. Pawlowicz 9/01 - Vectorized, changed functionality, put into standard M_MAP form.
%
% 29/Sep/2005 - fixed bug (apparently no-one ever used this option before now!)
%
% 1 Nov 2012 - added IGRF2011 coords (thanks to Joe Kinrade)
%
% Jan 2 2019 - Added IGRF2020 coords because motions have sped up a lot!
% Jan/2020 - Added IGRF option using user-specified date (thanks to W.
%            Brown)
% Jan/2020  - m_mag2geo and m_geo2mag did not work properly, they do now
%             (Added some more cases).

global MAP_PROJECTION MAP_COORDS

pi180=pi/180;

if nargin>=2 && isstruct(longIN)
  latIN=longIN.mdate;
  longIN=longIN.name;
end

switch optn
  case 'name'

    longOUT.name={'geographic','IGRF2000-geomagnetic','IGRF2011-geomagnetic','IGRF2020-geomagnetic','IGRF-geomagnetic'};
    return;
    
  case 'parameters'
    switch longIN
        case 'geographic'
            longOUT.name.name='geographic';
            longOUT.name.mdate=0;
            longOUT.lambda = 0;
            longOUT.phi = 0;
        
        case 'IGRF2000-geomagnetic'
            g10=-29615;
            g11=-1728;
            h11=5186;
            longOUT.name.name='IGRF2000-geomagnetic';
            longOUT.name.mdate=datenum(2000,1,1);
            longOUT.lambda=  atan(h11/g11);         
            longOUT.phi   =  atan((g11*cos(longOUT.lambda)+h11*sin(longOUT.lambda))/g10);

        case 'IGRF2011-geomagnetic' %ADDED 2012-11-02 Joe Kinrade, Uni. of Bath
            g10=-29496.5; %1st order IGFR2011 coefficient
            g11=-1585.9;  %2nd order
            h11=4945.1;   %3rd order
            longOUT.name.name='IGRF2011-geomagnetic';
            longOUT.name.mdate=datenum(2011,1,1);
            longOUT.lambda=  atan(h11/g11);         
            longOUT.phi   =  atan((g11*cos(longOUT.lambda)+h11*sin(longOUT.lambda))/g10);
            
        case 'IGRF2020-geomagnetic'  % ADDED 2020-01-03 W. Brown, British Geological Survey
           g10=-29404.8; % from IGRF-13 @ 2020/1/1 00:00:00
           g11=-1450.9;
           h11=4652.5;
           longOUT.name.name='IGRF2020-geomagnetic';
           longOUT.name.mdate=datenum(2020,1,1);
           longOUT.lambda=  atan(h11/g11);         
           longOUT.phi   =  atan((g11*cos(longOUT.lambda)+h11*sin(longOUT.lambda))/g10);
            
        case 'IGRF-geomagnetic'  %  Function that works for an arbitrary date.
                                 %  between 1900 and (predicted) to 2025 -
                                 %  see  mc_igrf.m for more details
           if nargin<3 || ~isfinite(latIN) 
               error(' You must explicitly specify a date when using the ''IGRF-geomagnetic'' option');
           end
           [g10,g11,h11]=mc_igrf(latIN);
           longOUT.name.name='IGRF-geomagnetic';
           longOUT.name.mdate=latIN;
           longOUT.lambda=  atan(h11/g11);         
           longOUT.phi   =  atan((g11*cos(longOUT.lambda)+h11*sin(longOUT.lambda))/g10);
            
      otherwise	
        error('Unrecognized coordinate system');
    end
    return;
    
  % These two cases are used in m_ll2xy and m_xy2ll ONLY
  case 'XYgeo2mag'   %  
     % Rotation matrices
     lambda=MAP_PROJECTION.coordsystem.lambda;
     phi=MAP_PROJECTION.coordsystem.phi;
     Tz=[ cos(lambda) sin(lambda) 0 ;
         -sin(lambda) cos(lambda) 0 ;
	      0       0           1 ];
		 
     Ty=[ cos(phi)  0   -sin(phi) ;
            0       1      0     ;
          sin(phi)  0   cos(phi) ]; 
     T=Ty*Tz;

  case 'XYmag2geo'
      % Rotation matrices
     lambda=MAP_COORDS.lambda;
     phi=MAP_COORDS.phi;
     Tz=[ cos(-lambda) sin(-lambda) 0 ;
         -sin(-lambda) cos(-lambda) 0 ;
	      0           0         1 ];
		 
     Ty=[ cos(-phi)  0   -sin(-phi) ;
            0        1         0     ;
          sin(-phi)  0    cos(-phi) ]; 
     T=Tz*Ty; % Reverse order
     
   % These two cases are used in m_mag2geo and m_geo2mag ONLY, and 
   % assume that an M_COORD call has ben made to set MAP_COORDS
      
    case 'geo2mag'
      % Rotation matrices
     lambda=MAP_COORDS.lambda;
     phi=MAP_COORDS.phi;
     Tz=[ cos(lambda) sin(lambda) 0 ;
         -sin(lambda) cos(lambda) 0 ;
	      0           0         1 ];
		 
     Ty=[ cos(phi)  0   -sin(phi) ;
            0        1         0     ;
          sin(phi)  0    cos(phi) ]; 
     T=Ty*Tz;

     case 'mag2geo'
      % Rotation matrices
     lambda=MAP_COORDS.lambda;
     phi=MAP_COORDS.phi;
     Tz=[ cos(-lambda) sin(-lambda) 0 ;
         -sin(-lambda) cos(-lambda) 0 ;
	      0           0         1 ];
		 
     Ty=[ cos(-phi)  0   -sin(-phi) ;
            0        1         0     ;
          sin(-phi)  0    cos(-phi) ]; 
     T=Tz*Ty; % reverse order

end
  
[n,m]=size(latIN);

% Degrees to radians, and make into rows.
latIN=(90 - latIN(:)')*pi180;
longIN=longIN(:)'*pi180;



% Transformation to cartesian coordinates.
x=sin(latIN).*cos(longIN);
y=sin(latIN).*sin(longIN);
z=cos(latIN);

% Rotations of the coordinates.
tmp=T*[x; y; z];
xp=tmp(1,:);yp=tmp(2,:);zp=tmp(3,:);

% Transformation back to spherical coordinates. Choosing the correct quadrant.
latOUT=acos(zp./sqrt(xp.^2+yp.^2+zp.^2));
longOUT=atan2(yp,xp);
	
if nargin > 3
	% Sign change due to sign convention used here.
	thetaVecIN=-thetaVecIN(:)';
        phiVecIN=phiVecIN(:)';
	
	% Computing vector rotations.
	% First, transformation to cartesian coordinates.
	
	% Rotation matrix is
	% [x] [  sLcG cLcG -sG ][radial]
	% [y]=[  sLsG cLsG  cG ][south ]
	% [z] [   cL   -sL  0  ][east  ]
	
	% Radial component is zero.
	Xvec= 0+    cos(latIN).*cos(longIN).*thetaVecIN - sin(longIN).*phiVecIN;
        Yvec= 0+    cos(latIN).*sin(longIN).*thetaVecIN + cos(longIN).*phiVecIN;
        Zvec= 0+                -sin(latIN).*thetaVecIN  + 0;
		
	% Rotations of the system.
	tmp=T*[Xvec;Yvec; Zvec];
	Xvecp=tmp(1,:);Yvecp=tmp(2,:);Zvecp=tmp(3,:);
	
	% Transformation back to spherical coordinates.

	thetaVecOUT=cos(latOUT).*cos(longOUT).*Xvecp + cos(latOUT).*sin(longOUT).*Yvecp -sin(latOUT).*Zvecp ;
	phiVecOUT  =            -sin(longOUT).*Xvecp +              cos(longOUT).*Yvecp              +0       ;
	
	% Sign change due to sign convention used here.
	thetaVecOUT=-reshape(thetaVecOUT,n,m);
        phiVecOUT  = reshape(phiVecOUT,n,m);
end


% Radians to degrees.
latOUT=90 - reshape(latOUT,n,m)/pi180;
longOUT=reshape(longOUT,n,m)/pi180;

