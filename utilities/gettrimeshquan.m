function tq = gettrimeshquan( p, t)
% tq.ds  - edge length connecting node
%           2-3, 3-1, and 1-2 respectively
% tq.J   - 0.5*element area
% tq.vang - angle assocaited with node 1, 2, 3 respectively
%                       
              
edge = [2  3  1
        3  1  2] ;

ne = size( t, 1) ;

p1 = p(t(:, edge(1,:))',:)  ;
p2 = p(t(:, edge(2,:))',:)  ;

%
ds = reshape(sqrt(sum((p1 - p2).^2,2)), 3, ne )' ;

idxcm = [ 1 2 3
          2 3 1
          3 1 2 ] ;

% Internal angle 
tq.vang = zeros(ne,3) ;
for i = 1: 3
   tq.vang(:,i) = acos((-ds(:,idxcm(i,1)).^2 +  ds(:,idxcm(i,2)).^2 + ...
           ds(:,idxcm(i,3)).^2)./(2*ds(:,idxcm(i,2)).*ds(:,idxcm(i,3)))) ;
end
%

% Length 
tq.ds = ds ;

% Jacobian
xr = 0.5*(p(t(:,2),:) - p(t(:,1),:)) ;
xs = 0.5*(p(t(:,3),:) - p(t(:,1),:)) ; 

tq.J = (xr(:,1).*xs(:,2) - xs(:,1).*xr(:,2)) ;

% idxtb = find( vang' < 25*pi/180 ) ;
% ietb = unique(ceil(idxtb/3)) ;

%
%  Bank, Randolph E., 
%      PLTMG: A Software Package for Solving Elliptic Partial Differential Equations, 
%      User's Guide 6.0, Society for Industrial and Applied Mathematics, Philadelphia, PA, 1990.
tq.qm = 4*sqrt(3)*(2*tq.J)./(sum(tq.ds.^2,2)) ; 