function fem = update3nbr(fem_struct,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This routine will remove an interior node that has only three nodal
% neighbors.  Along with the node being removed, two elements will be
% removed and the remaining element will be updated. The node will be
% removed from the fem_struct.x,y,z,nei. The elements will be removed
% from the fem_struct.e and fem_struct.ar will be corrected.
%
% Additionally, fem_struct.e, fem_struct.nei, and fem_struct.bnd
% will be updated to that its components point to the new
% numbering of the nodes.
%
%            OLD                   NEW
%              +                     +
%             /|\                   / \
%            / | \                 /   \
%           /  |  \               /     \
%          /   *   \             /       \
%         /   / \   \           /         \
%        /  /     \  \         /           \
%       / /         \ \       /             \
%      +---------------+     +---------------+
%
%
%  Use this one with caution. 
%  I would recommend only using it after all the 
%  update 8,9,etc... have been updated,
%  and in such a way that only nodes that weren't 4 before
%  be updated.  Also, don't use it on boundary nodes.
%
%
%  Usage:
%        fem_struct = update3nbr(fem_struct,i)
%
%        fem_struct -- the OPNML finite element structure
%        i -- the node number to be eliminated
%
%
%  Name:        update3nbr.m
%  Written By:  Chris Massey
%  Date:        Feb. 2, 2005 
%  Last Modifed:  
%         Sept. 19, 2006 -- Chris Massey, NRL Code 7322, S.S.C. MS 39529
%                      % Add test to make sure the correct nodal neighbor
%                        connectivity existed for the requested update
%                        and added a test to create the nodal neighbor list
%                        if not present.
%         Nov. 30, 2007 -- Chris Massey, NRL Code 7322, SSC, MS
%                    %Fixed indexing bug, changed jj to i in line 56-60,
%                    % Fixed bug in which nodes to keep in the list and
%                      added updating of the boundary node list.
%         July 14, 2008 -- Chris Massey, NRL Code 7322
%               moved test for 4 neighbors to after the test for nei.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keywords:  reduce connectivity, update 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checks if the fem_struct has an nei. If not one is generated.
try
    nei = fem_struct.nei;
catch
    fem_struct.nei = ele2nei(fem_struct.e,fem_struct.x,fem_struct.y);
    nei = fem_struct.nei;
end

tempnc = sum(~ismember((fem_struct.nei(i,:)),0));
if tempnc ~= 3
    error(['There are ',num2str(tempnc),' nodal neighbors for node number ',...
           num2str(i),'. This routine only works for 3 nodal neighbors.']);
end

fem = fem_struct;
x = fem_struct.x;
y = fem_struct.y;
z = fem_struct.z;
enodes = fem_struct.e;
elem_area = fem_struct.ar;
bnd = fem_struct.bnd;

angmat = [4,1,2;1,2,3;2,3,4;3,4,1];
vod = [1,2,3,4,1,2,3,4];


% Find all the element numbers that contain the node to 
% be destroyed
[I2,J2] = find(enodes == i);
clear J2

% Elements to keep & to eliminate
elem_keep = I2(1);
elem_toss = I2([2,3]);

% Update the area for the elements being kept
elem_area(elem_keep) = sum(elem_area(I2));

% Sort the tossed elements for later use
elem_toss = sort(elem_toss);

% Find the nodes that define the outer big triangle 
%node_perim(1) = setdiff(enodes(elem_keep,:),enodes(elem_toss(1),:)); 
%node_perim(2:3) = setxor(enodes(elem_keep,:),i); % The other nodes in the 
%                                            % element being kept
node_perim = nei(i,1:3);                                            
                                            
% Update the element that is being kept
% by replacing the destroyed node with the 
% correct outer node
I3 = find(enodes(elem_keep,:)==i);
enodes(elem_keep,I3) = setdiff(node_perim,enodes(elem_keep,:));

%Update Element List (Destory two elements)
enodes = enodes([1:elem_toss(1)-1,elem_toss(1)+1:elem_toss(2)-1,elem_toss(2)+1:end],:);
elem_area = elem_area([1:elem_toss(1)-1,elem_toss(1)+1:elem_toss(2)-1,elem_toss(2)+1:end]);

%Update nei list
for j=1:3
    I3 = find(nei(node_perim(j),:)==i);
    nei(node_perim(j),[I3:end]) = [nei(node_perim(j),[I3+1:end]),0];
end


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

I3 = find(bnd > i);
bnd(I3) = bnd(I3)-1;

% Define the new fem structure
fem = fem_struct;
fem.x = x;
fem.y = y;
fem.z = z;
fem.nei = nei;
fem.e = enodes;
fem.bnd = bnd;
fem.ar = elem_area;

return
