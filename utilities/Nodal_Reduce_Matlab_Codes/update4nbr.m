function fem = update4nbr(fem_struct,i)
% update4nbr removes an interior node that has
%  only four nodal neighbors.  Along with the node
%  being removed, two elements will be removed and 
%  the remaining two elements will be updated.
%  The node will be removed from the fem_struct.x,y,z,nei
%  The elements will be removed from the fem_struct.e
%  and fem_struct.ar will be corrected.
%
% Additionally, fem_struct.e will be updated to that 
%  its components point to the new numbering of the nodes.
%
%    (Original)             (Opt.A)       (Opt.B)
%     X-----X               X-----X       X-----X
%     |\   /|               |\    |       |    /|
%     | \ / |               | \   |       |   / |
%     |  X  |   becomes     |  \  |  or   |  /  |
%     | / \ |               |   \ |       | /   |
%     |/   \|               |    \|       |/    |
%     X-----X               X-----X       X-----X
%
% Opt. A or Opt. B is determined by looking at all the
% surrounding nodes and determining the one with the 
% largest area.  The one with the largest area is selected
% to have it's diagonal to the old node left in place.
%
%
%  Use this one with caution. 
%  I would recommend only using it after all the 
%  update 8,9,etc... have been updated,
%  and in such a way that only nodes that weren't 4 before
%  get updated.  Also, don't use it on boundary nodes.
%
%  Usage:
%        fem_struct = update4nbr(fem_struct,i)
%
%        fem_struct -- the OPNML finite element structure
%        i -- the node number to be eliminated
%
%
%  Name:        update4nbr.m
%  Written By:  Chris Massey
%  Date:        May 6, 2004
% Last Modifed:  
%         Sept. 19, 2006 -- Chris Massey, NRL Code 7322, S.S.C. MS 39529
%                      % Add test to make sure the correct nodal neighbor
%                        connectivity existed for the requested update
%                        and added a test to create the nodal neighbor list
%                        if not present.
%         April 11, 2006 -- Chris Massey, NRL Code 7322
%                   Fixed a typo bug.
%         July 14, 2008 -- Chris Massey, NRL Code 7322
%               moved test for 4 neighbors to after the test for nei.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Checks if the fem_struct has an nei. If not one is generated.
try
    nei = fem_struct.nei;
catch
    fem_struct.nei = ele2nei(fem_struct.e,fem_struct.x,fem_struct.y);
    nei = fem_struct.nei;
end

tempnc = sum(~ismember((fem_struct.nei(i,:)),0));
if tempnc ~= 4
    error(['There are ',num2str(tempnc),' nodal neighbors for node number ',...
           num2str(i),'. This routine only works for 4 nodal neighbors.']);
end

fem = fem_struct;
x = fem_struct.x;
y = fem_struct.y;
z = fem_struct.z;
enodes = fem_struct.e;
elem_area = fem_struct.ar;

angmat = [4,1,2;1,2,3;2,3,4;3,4,1];
vod = [1,2,3,4,1,2,3,4];


% Compute the Quadrilateral's Angles
for j = 1:4
    a = sqrt((x(nei(i,angmat(j,1)))-x(nei(i,angmat(j,2))))^2+ (y(nei(i,angmat(j,1)))-y(nei(i,angmat(j,2))))^2);
    b = sqrt((x(nei(i,angmat(j,2)))-x(nei(i,angmat(j,3))))^2+ (y(nei(i,angmat(j,2)))-y(nei(i,angmat(j,3))))^2);
    c2 = (x(nei(i,angmat(j,1)))-x(nei(i,angmat(j,3))))^2+(y(nei(i,angmat(j,1)))-y(nei(i,angmat(j,3))))^2;
    quadang(j) = acos((a^2+b^2-c2)/(2*a*b));
end
%quadang = quadang*180/pi

% Find the largest angle
I = find(quadang == max(quadang));
I = I(1);

% Find all the element numbers that contain the node to 
% be destroyed
[I2,J2] = find(enodes == i);
clear J2

% Local copy of the "enodes" for the elements involved
for j = 1:4
    temp = ismember(sort(enodes(I2,:),2),sort([i,nei(i,vod(j)),nei(i,vod(j+1))]),'rows');
    I3 = find(temp == 1);
    elem(j,1:3) = enodes(I2(I3),1:3);
end

% Elements to keep & to eliminate
[elem_keep,J2] = find(enodes(I2,:) == nei(i,I));
elem_keep = I2(elem_keep);
elem_toss = setdiff(I2,elem_keep);
%elem_keep = I2([vod(I),vod(I+1)]);
%elem_toss = I2([vod(I+2),vod(I+3)]);


% Update the area for the elements being kepted
elem_area(elem_keep) = elem_area(elem_keep)+elem_area(elem_toss);

% Sort the tossed elements for later use
elem_toss = sort(elem_toss);

% Find the off diagonal nodes and the other node
% that stays on the diagonal
off_diag_nodes = nei(i,angmat(I,[1,3]));
node_update = setdiff(nei(i,1:4),[i,nei(i,angmat(I,1:3))]);
diag_nodes = setdiff(unique(enodes([elem_keep,elem_toss],:)),[i,off_diag_nodes]);

% Update the elements that are being kepted
% by replacing the destroyed node with the 
% new node
for j = 1:2
    I3= find(enodes(elem_keep(j),:)==i);
    %enodes(elem_keep(j),:)
    enodes(elem_keep(j),I3) = node_update;
    %enodes(elem_keep(j),:)
end


%Update Element List
enodes = enodes([1:elem_toss(1)-1,elem_toss(1)+1:elem_toss(2)-1,elem_toss(2)+1:end],:);
elem_area = elem_area([1:elem_toss(1)-1,elem_toss(1)+1:elem_toss(2)-1,elem_toss(2)+1:end]);

%Update nei list
for j = 1:2
    I3 = find(nei(off_diag_nodes(j),:)==i);
    nei(off_diag_nodes(j),[I3:end]) = [nei(off_diag_nodes(j),[I3+1:end]),0];
end
I3 = find(nei(diag_nodes(1),:) == i);
nei(diag_nodes(1),I3) = diag_nodes(2);
I3 = find(nei(diag_nodes(2),:) == i);
nei(diag_nodes(2),I3) = diag_nodes(1);
    

% Delete the components associated with the destroyed 
% node
nei = nei([1:i-1,i+1:end],:); %Deletes the i(th) nei
x = x([1:i-1,i+1:end]);
y = y([1:i-1,i+1:end]);
z = z([1:i-1,i+1:end]);

% Update the elements list
% by using the new node numbers
I3 = find(enodes > i);
enodes(I3) = enodes(I3)-1;

I3 = find(nei > i);
nei(I3) = nei(I3)-1;

I3 = find(fem_struct.bnd > i);
fem_struct.bnd(I3) = fem_struct.bnd(I3)-1;

% Define the new fem structure
fem = fem_struct;
fem.x = x;
fem.y = y;
fem.z = z;
fem.e = enodes;
fem.ar = elem_area;
fem.nei = nei;

