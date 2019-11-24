function [pc,tc]=collapse_thin_triangles(p,t,minqm)
% Collapse triangles that are exceedingly thin incrementally modifying the 
% connectivity.
% Thin triangles are identified by containing highly actue angles.
% The shortest edge of the identified triangle is collapsed to a point.
% and the triangle table and point matrix are updated accordingly.
% kjr, usp,br. 2019.

if nargin < 3
    minqm = 0.10 ; % minimum quality to identify thin triangles
end
tq = gettrimeshquan(p,t);

kount=1;
while any(tq.qm < minqm)
    
    id = find(tq.qm < minqm); % queue

    tid=id(1);
    tlocal=t(tid,:);
    ee = [tlocal(1,[1,2]); tlocal(1,[1,3]); tlocal(1,[2,3])];
    evec = p(ee(:,1),:)- p(ee(:,2),:);
    [~,meid]=min(sum(evec.^2,2)); % find shortest edge
    
    PointIDToRemove=ee(meid,1) ; % remove this id for ReplaceID 
    PointIDToReplace=ee(meid,2); % replace remove ID with this id 
 
    % Delete triangle that has the thin edge for deletion
    t(tid,:)=[]; 
    % Replace all other instances of PointIDToRemove with PointIDToReplace 
    t(t==PointIDToRemove)=PointIDToReplace;  

     % recompute qualities. 
    tq = gettrimeshquan(p,t);
    kount=kount+1; 
end
[pc,tc]=fixmesh(p,t); % remove hanging points 
disp(['Removed ',num2str(kount),' thin triangles!']); 


