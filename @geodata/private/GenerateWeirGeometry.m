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

tmp = crestline ; crestline = [] ; 

[crestline(:,2),crestline(:,1)] = my_interpm(tmp(:,2),tmp(:,1),spacing) ; 

width = width/2 ; 

dy = gradient(crestline(:,2));
dx = gradient(crestline(:,1));

% % plt the tangent vec
% quiver(crestline(:,1),crestline(:,2),-dy,dx)
% hold on; plot( crestline(:,1), crestline(:,2))

% tangential vectors to line points
tang(:,1) = crestline(:,1)-dy ;
tang(:,2) = crestline(:,2)+dx ;

v=tang-crestline;

u=v./norm(v) ;

new_above=crestline + width*u ;
new_below=crestline - width*u ;

% create ibconn_pts 
ibconn_pts = [new_above new_below] ; 
% hold on; plot(new_above(:,1),new_above(:,2),'r-s') ;
% hold on; plot(new_below(:,1),new_below(:,2),'g-s') ;

% create knife-edge corners by calculating the normal of the new sgements
% formed at the edges of the geometry
ke_left = [new_above(1,:); new_below(1,:)] ;

dy = gradient(ke_left(:,2));
dx = gradient(ke_left(:,1));

% tangential vectors to line points
tang2(:,1) = ke_left(:,1)-dy ;
tang2(:,2) = ke_left(:,2)+dx ;

v=tang2-ke_left;

u=v./norm(v) ;

ke_left_point=ke_left - 4.0*width*u ;
ke_left_point_mid = mean(ke_left_point) ;

%hold on; plot(ke_left_point_mid(:,1),ke_left_point_mid(:,2),'b*') ;

% create knife-edge corners by calculating the normal of the new sgements
% formed at the edges of the geometry
ke_right = [new_above(end,:); new_below(end,:)] ;

dy = gradient(ke_right(:,2));
dx = gradient(ke_right(:,1));

% tangential vectors to line points
tang3(:,1) = ke_right(:,1)-dy ;
tang3(:,2) = ke_right(:,2)+dx ;

v=tang3-ke_right;

u=v./norm(v) ;

ke_right_point=ke_right + 4.0*width*u ;
ke_right_point_mid = mean(ke_right_point) ;

%hold on; plot(ke_right_point_mid(:,1),ke_right_point_mid(:,2),'b*') ;

weirPfix = [ke_left_point_mid ; new_below; ke_right_point_mid ;  flipud(new_above) ; ke_left_point_mid] ;

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



