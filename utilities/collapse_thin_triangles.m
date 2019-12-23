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
tqm = tq.qm;
kount = 0;
while any(tqm < minqm)
    
    id = find(tqm < minqm); % queue

    tid=id(1);
    tlocal=t(tid,:);
    ee = [tlocal(1,[1,2]); tlocal(1,[1,3]); tlocal(1,[2,3])];
    evec = p(ee(:,1),:)- p(ee(:,2),:);
    [~,meid]= min(sum(evec.^2,2)); % find shortest edge
    
    PointIDToRemove=ee(meid,1) ; % remove this id for ReplaceID 
    PointIDToReplace=ee(meid,2); % replace remove ID with this id 
    
    % Delete triangle that has the thin edge for deletion
    to = t;
    t(tid,:)=[]; 
    % Replace all other instances of PointIDToRemove with PointIDToReplace
    t(t==PointIDToRemove) = PointIDToReplace;  
    I = t(:,1) - t(:,3) == 0 | t(:,2) - t(:,3) == 0 | t(:,1) - t(:,2) == 0; 
    if sum(I) > 0
        % cannot delete this element, make the quality good.
        t = to; tqm(tid) = 1;
        continue;
    end
     % recompute qualities. 
    tqm(tid) = [];
    I = t(:,1) == PointIDToReplace | t(:,2) == PointIDToReplace | t(:,3) == PointIDToReplace;
    tq = gettrimeshquan(p,t(I,:));
    tqm(I) = tq.qm;
    kount=kount+1; 
end
[pc,tc]=fixmesh(p,t); % remove hanging points 
disp(['Removed ',num2str(kount),' thin triangles!']); 


