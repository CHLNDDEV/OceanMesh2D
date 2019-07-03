function [weirPfix,weirEgfix,ibconn_pts] = GenerateWeirGeometry(crestline,width,spacing,plot_on)
% Builds weir(d) island geometries from a weir
% crestline polyline and returns the points in a order that define the
% boundary.
%
% INPUTS:
%
% crestline: a vector of points that define the crest of the weir
%
% width: a scalar that defins the width (in WGS84 degrees)
%
% spacing: ensure point spacing is bounded above by spacings  
%
% plot_on: 1 for option to visualize faux weir island geometry
%
% OUTPUTS:
%
% weirPoints: points of the weir faux island boundary
%
% weirEdges:  edges that index into weirPoints defining the connections
%             between points on the weir island boundary
%
% ibconn_pts: This array tells you the coordinates that form a pair of points.  
%
% AUTHOR: KJR, OCT 26. 2018, CHL,UND
%    updated by KJR May, 23, 2019, CHL,UND first and last point correspond to knife
%    edges that are the minimum element size long 


tmp = crestline ; crestline = [] ; 

[crestline(:,2),crestline(:,1)] = my_interpm(tmp(:,2),tmp(:,1),spacing) ; 

width = width/2 ; 

dy = gradient(crestline(:,2));
dx = gradient(crestline(:,1));

% tangential vectors to line points
tang(:,1) = crestline(:,1)-dy ;
tang(:,2) = crestline(:,2)+dx ;

% normal
v=tang-crestline;
% unit normal 
u=v./norm(v) ;

new_above=crestline(2:end-1,:) + width*u(2:end-1,:) ;
new_below=crestline(2:end-1,:) - width*u(2:end-1,:) ;

% create ibconn_pts 
ibconn_pts = [new_above new_below] ; 
ke_right_point_mid = crestline(end,:) ; 
ke_left_point_mid  = crestline(1,:) ; 

weirPfix = [ke_left_point_mid ; new_below; 
            ke_right_point_mid ; 
           flipud(new_above) ; ke_left_point_mid] ;

% compute the edges that define the weir geometry for output
weirEgfix = Get_poly_edges([weirPfix; NaN NaN]) ;
weirPfix(end,:) = []; % <--this point was dup'ped
change=find(weirEgfix==max(weirEgfix));
weirEgfix(change) = weirEgfix(change)-1 ; 
if nargin > 2
    if plot_on 
        figure;
        %simpplot(weirIsland,tweir) ;
        hold on ; plot(weirPfix(:,1),weirPfix(:,2),'r-');
        hold on; plot(crestline(:,1),crestline(:,2),'kx-','linewi',2) ;
        axis off
    end
end
end



