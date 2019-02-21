function [Fx,Fy] = comptide_eflux( b, eamp, ephs, uamp, uphs, vamp, vphs, bcut )
% calculate the tidal enegry flux averaged 
% over one tidal cycle:
%
%   (Fx,Fy) = b (1/T) \int_{0}^{T} (u \eta, v \eta) dx
%
%    b = bathymetry
%    eta = eamp*cos( wt - ephs )
%    u = uamp*cos( wt - uphs )
%    v = vamp*cos( wt - vphs )
%
% Note:
%  Amplitude in meter
%  Phase in  radian
%
% Return:  Fx/(\rho*g) 
%          Fy/(\rho*g)
Fx = 0.5*b.*(uamp.*eamp).*cos( uphs - ephs ) ; 
Fy = 0.5*b.*(vamp.*eamp).*cos( vphs - ephs ) ; 

bc = 0 ; 
if ( nargin > 7 )
    bc = bcut ;
end
idx = find( b < bc ) ; % Land 
if ( ~isempty( idx ) ) 
  Fx(idx) = 0.0 ;
  Fy(idx) = 0.0 ;
end

