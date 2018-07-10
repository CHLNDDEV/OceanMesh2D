function fem = patch_update(oldfem,newfem,jj);
% PATCH_UPDATE - retriangulates the newly placed points from the
% update#nbr routines and incorporates that triangulation into
% the new fem struct.
%
% Variable
%  jj -- the node number from the old mesh that is being updated
%  oldfem -- the old fem grid structure with high connectivity
%  newfem -- the new fem grid structure with lowered connectivity
%  fem -- a new triangulated fem with hopefully lower connectivity
%
%  Calls:
%  orientccw.m, ele2nei.m
%
%
%  FileName:    patch_update.m
%  Written By:  Chris Massey
%  Date:        Aug.  5, 2005
%  Modified:    Nov.  1, 2005 -- Fixed bug
%               Apr. 17, 2007 -- Chris Massey, NRL Code 7322
%                 Added section to ensure all the new created points
%                 and modified original points are still in the patch
%                 before attempting to triangulate.
%               April. 23, 2007 -- Chris Massey,
%                 Added section to try triangulating original patch 
%                 configuration if the new configuration failed.  If both
%                 new triangulations fail, then we now return the original
%                 mesh.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nnhbr = find(oldfem.nei(jj,:)==0,1,'first');
if isempty(nnhbr) == 0
    nnhbr = nnhbr - 1;
else
    nnhbr = size(oldfem.nei,2);
end


nnold = length(oldfem.x);
nnnew = length(newfem.x);
neold = size(oldfem.e,1);
nenew = size(newfem.e,1);

node_list = [oldfem.nei(jj,1:nnhbr),jj,nnold+1:nnnew];
[I,J] = find(oldfem.e == jj);
elem_list = sort([I',neold+1:nenew]);
nepatchnew = length(elem_list);

clear I J


x0 = newfem.x(jj);y0 = newfem.y(jj);z0 = newfem.z(jj);
xex = oldfem.x(oldfem.nei(jj,1:nnhbr));
yex = oldfem.y(oldfem.nei(jj,1:nnhbr));
zex = oldfem.z(oldfem.nei(jj,1:nnhbr));
xin = newfem.x(nnold+1:nnnew);
yin = newfem.y(nnold+1:nnnew);
zin = newfem.z(nnold+1:nnnew);

% TCM 04/17/2007 -- Begin Changes
% Make sure all the new points (and the old center point)
% lie inside the patch.  If not move each outlying point
% closer to the original center location (midpoint between
% original center and current location).
xt1 = [x0;xin(:)];
yt1 = [y0;yin(:)];
zt1 = [z0;zin(:)];
[in,on] = inpolygon(xt1,yt1,xex,yex);
I = find(in == 0 | on == 1);
while ~isempty(I)
    xt1(I) = (xt1(I) + oldfem.x(jj))./2;
    yt1(I) = (yt1(I) + oldfem.y(jj))./2;
    zt1(I) = (zt1(I) + oldfem.z(jj))./2;
    [in,on] = inpolygon(xt1,yt1,xex,yex);
    I = find(in == 0 | on == 1);
end
x0 = xt1(1);y0 = yt1(1);z0 = zt1(1);
xin = xt1(2:end);yin = yt1(2:end);zin = zt1(2:end);
clear xt1 yt1 zt1 in on
%
% TCM 04/17/2007 -- End Changes

patchfem = patch_tri(x0,y0,z0,xex,yex,zex,xin,yin,zin);
nepatch = size(patchfem.e,1);

if(nepatch ~= 0) % Newly triangulated patch with new nodes is valid
    fem = newfem;
    fem.x(node_list) = patchfem.x;
    fem.y(node_list) = patchfem.y;
    fem.z(node_list) = patchfem.z;

    tempfem = patchfem;
    tempfem = orientccw(tempfem);

    for i=1:length(node_list)
        I = find(patchfem.e == i);
        tempfem.e(I) = node_list(i);
        clear I
    end


    patchfem = tempfem;
    clear tempfem;

    %Update the fem.e with the new connectivity list

    %Number of elements new patch is smaller than or equal to old patch number
    % of elements
    if nepatch<=nepatchnew
        for i=1:nepatch
            fem.e(elem_list(i),:) = patchfem.e(i,:);
        end
        if nepatch < nepatchnew
            keepelist = setdiff(1:1:nenew,elem_list(nepatch+1:end));
            fem.e = fem.e(keepelist,1:3);
            clear keepelist
        end
    else
        % Number of elements in new patch is larger than the number in the
        % old patch
        for i=1:nepatchnew
            fem.e(elem_list(i),1:3) = patchfem.e(i,1:3);
        end
        pb = 1;
        for i=nepatchnew+1:nepatch
            fem.e(nenew+pb,1:3) = patchfem.e(i,1:3);
            pb = pb + 1;
        end
    end

    fem = orientccw(fem);
    fem.nei = ele2nei(fem.e,fem.x,fem.y);

else % Newly triangulated patch with new nodes is invalid,
     % Try triangulating old patch and see if it is valid.
     
    patchfem = patch_tri(oldfem.x(jj),oldfem.y(jj),oldfem.z(jj),...
        xex,yex,zex,[],[],[]);
    nepatch = size(patchfem.e,1);
    if nepatch ~=0 % The newly triangluated patch is valid
        tempfem = patchfem;
        tempfem = orientccw(tempfem);

        node_list = [oldfem.nei(jj,1:nnhbr),jj];
        fem = oldfem;
        fem.x(node_list) = patchfem.x;
        fem.y(node_list) = patchfem.y;
        fem.z(node_list) = patchfem.z;

        [I,J] = find(oldfem.e == jj);
        elem_list = sort([I']);
        nepatchnew = length(elem_list);
      
        for i=1:length(node_list)
            I = find(patchfem.e == i);
            tempfem.e(I) = node_list(i);
            clear I
        end
        patchfem = tempfem;
        clear tempfem;
        %Number of elements new patch is smaller than or equal to old patch number
        % of elements
        if nepatch<=nepatchnew
            
            for i=1:nepatch
                fem.e(elem_list(i),:) = patchfem.e(i,:);
            end
            if nepatch < nepatchnew
                keepelist = setdiff(1:1:neold,elem_list(nepatch+1:end));
                fem.e = fem.e(keepelist,1:3);
                clear keepelist
            end
        else
            % Number of elements in new patch is larger than the number in the
            % old patch
            for i=1:nepatchnew
                fem.e(elem_list(i),1:3) = patchfem.e(i,1:3);
            end
            pb = 1;
            for i=nepatchnew+1:nepatch
                fem.e(neold+pb,1:3) = patchfem.e(i,1:3);
                pb = pb + 1;
            end
        end

        fem = orientccw(fem);
        fem.nei = ele2nei(fem.e,fem.x,fem.y);


    else % Retriangulating failed in both cases (new nodes, old nodes)

        fem = oldfem;  % The newly triangulated mesh did not conform to the
        % original patch boundary.
        disp(' ');
        disp('The newly triangulated mesh did not conform to the');
        disp('original patch boundary.  Returning original mesh.');
        disp(' ');
    end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Private function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fem] = patch_tri(x0,y0,z0,xex,yex,zex,xin,yin,zin)
% PATCH_TRI takes a set of exterior and interior points on and in a patch
% and calls delaunay triganulation then creates a small fem struct of
% that newly triangluated patch.
%
% Variables
% x0,y0,z0 -- the coordinates of the center of the original patch
% xex,yex,zex -- the coordinates of the exterior of the original patch
% xin,yin,zin -- the coordinates of the added interior points of the patch
%
% Calls:
%  delaunay.m, detbndy.m, springs.m, orientccw.m
%
%  Name:   patch_tri.m
%  Written By:  Chris Massey
%  Date:        Aug. 5, 2005
% Modified:     Oct. 18, 2005 -- fixed possible logic flaw in
%                                removing newly created boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [xex;x0;xin];
y = [yex;y0;yin];
z = [zex;z0;zin];


% Triangulate the points and create a fem_grid structure
tri = delaunay(x(:),y(:));%  ,{'Qt','Qbb','Qc','Qz'}); -- removed qhull specific options
fem_struct.name='patch';
fem_struct.x = x(:);
fem_struct.y = y(:);
fem_struct.z = z(:);
fem_struct.e = tri;
fem_struct.bnd=detbndy(fem_struct.e);

% Create a boundary list for the original patch based on the
% exterior nodes.  Sort that boundary list and the newly
% trianglulated boundary list

bnd1(:,1) = (1:1:length(xex))';
bnd1(:,2) = [2:1:length(xex),1];
bnd1 = sort(bnd1,2);
bnd2 = sort(fem_struct.bnd,2);

% Identify any newly created boundary segments
% and delete the corresponding elements
temp = setdiff(bnd2,bnd1,'rows');

while (isempty(temp) == 0 )

    for i=1:size(temp,1);
        [I3,J3] = find(fem_struct.e == temp(i,1));
        [I4,J4] = find(fem_struct.e == temp(i,2));
        I5 = intersect(I3,I4);
        fem_struct.e = fem_struct.e([1:I5-1,I5+1:end],:);
    end
    if ~isempty(fem_struct.e)
        fem_struct.bnd = detbndy(fem_struct.e);
        bnd2 = sort(fem_struct.bnd,2);
        temp = setdiff(bnd2,bnd1,'rows');
    else
        fem = fem_struct;
        temp = [];
        %Empty element list;
    end

end


% Make sure all nodes are listed in CCW,
% then spring the new region creating the final new structure.
if ~isempty(fem_struct.e)
    fem = fem_struct;
    fem = orientccw(fem_struct);
    fem = springs(fem);
end


return