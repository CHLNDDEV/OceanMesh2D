function fem = update11nbr(fem_struct,jj)
% update11nbr takes a node with 11 elements connected to it 
%  (which defines a region) and addes 4 nodes and 8 elements
%  to that region then regroups so that no node has more than
%  six elements connected to it. The new nodes and elements are sprung.
%
% NOTE:  fem_struct.nei must be a component.
%
%  Variables
%  fem_struct -- finite element structure from opnml.
%  jj -- the node number that has 11 elements connected to it.
%  fem -- the updated finite element structure.
%        Note:  Only fem.x,fem.y,fem.z,fem.e, and fem.nei
%               get updated.  fem.bnd does not need to be updated.
%
%  Usage -- fem = update11nbr(fem_struct,jj)
%
%  Name: update11nbr.m
%  Written by: Ben Holladay
%  Date: June 22,2004
%  Updated:  Aug. 27, 2004 -- Chris Massey Fixed bug in .nei for the
%                                          3rd new node.
% Last Modifed:  
%         Sept. 19, 2006 -- Chris Massey, NRL Code 7322, S.S.C. MS 39529
%                      % Add test to make sure the correct nodal neighbor
%                        connectivity existed for the requested update
%                        and added a test to create the nodal neighbor list
%                        if not present.
%         Aug. 13, 2007 -- Chris Massey, NRL Code 7322
%                     Fixed a bug in the code that updates nodal neighbor
%                     list for nodes getting extra connectivity.
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

tempnc = sum(~ismember((fem_struct.nei(jj,:)),0));
if tempnc ~= 11
    error(['There are ',num2str(tempnc),' nodal neighbors for node number ',...
           num2str(jj),'. This routine only works for 11 nodal neighbors.']);
end

%Sets the intial fem_struct variables.
x = fem_struct.x;
y = fem_struct.y;
z = fem_struct.z;
enodes = fem_struct.e;
vod = [1:11,1:11];

%Adds new nodes and elements and labels the bad node, neighboring nodes,
%and neighboring elements.
[nnodes,nc] = size(nei);
nelems = size(enodes,1);
newnode = nnodes + [1:4];
newelem = nelems + [1:8];
badnode = jj;
nbrnode = nei(jj,1:11);
[nbrelem,J] = find(enodes == badnode);
clear J;
inrad = max(sqrt((x(nbrnode)-x(badnode)).^2 + (y(nbrnode)-y(badnode)).^2));

%First guess of the new coordinates and bathymetry.
Bx = sum(x([(nbrnode([1:3])),badnode]))/4;
By = sum(y([(nbrnode([1:3])),badnode]))/4;
Bz = sum(z([(nbrnode([1:3])),badnode]))/4;

Ax = sum([(x([(nbrnode([10,11,1]))]))',Bx])/4;
Ay = sum([(y([(nbrnode([10,11,1]))]))',By])/4;
Az = sum([(z([(nbrnode([10,11,1]))]))',Bz])/4;

Cx = sum([(x([(nbrnode([3:6]))]))',Bx])/5;
Cy = sum([(y([(nbrnode([3:6]))]))',By])/5;
Cz = sum([(z([(nbrnode([3:6]))]))',Bz])/5;

Ex = sum(x([(nbrnode([7:10])),badnode]))/5;
Ey = sum(y([(nbrnode([7:10])),badnode]))/5;
Ez = sum(z([(nbrnode([7:10])),badnode]))/5;

Dx = sum([(x([(nbrnode([6:7]))]))',Bx,Ax,Cx,Ex])/6;
Dy = sum([(y([(nbrnode([6:7]))]))',By,Ay,Cy,Ey])/6;
Dz = sum([(z([(nbrnode([6:7]))]))',Bz,Az,Cz,Ez])/6;

x(badnode) = Ax;
y(badnode) = Ay;
z(badnode) = Az;
x(newnode(1)) = Bx;
y(newnode(1)) = By;
z(newnode(1)) = Bz;
x(newnode(2)) = Cx;
y(newnode(2)) = Cy;
z(newnode(2)) = Cz;
x(newnode(3)) = Dx;
y(newnode(3)) = Dy;
z(newnode(3)) = Dz;
x(newnode(4)) = Ex;
y(newnode(4)) = Ey;
z(newnode(4)) = Ez;

%Updates the effected elements. J variables are used to represent which
%elements are connected each node.
i = 1;
while i <= 9;
j1 = ismember((enodes((nbrelem),:)),(nbrnode(vod(i))));
j2 = ismember((enodes((nbrelem),:)),(nbrnode(vod(i+1))));
spb = [1,1,2,2,2,3,4,4,4];
j = sum([j1,j2],2);
temp = nbrelem(find (j == 2));
enodes(temp,:) = ([nbrnode(vod(i)),nbrnode(vod(i+1)),newnode(spb(i))]);
i = i + 1;
clear j j1 j2 temp;
end

%Addes the 8 new elements.
enodes(newelem(1),:) = [nbrnode(1),newnode(1),badnode];
enodes(newelem(2),:) = [nbrnode(3),newnode(2),newnode(1)];
enodes(newelem(3),:) = [newnode(1),newnode(2),newnode(3)];
enodes(newelem(4),:) = [newnode(1),newnode(3),badnode];
enodes(newelem(5),:) = [badnode,newnode(3),newnode(4)];
enodes(newelem(6),:) = [nbrnode(10),badnode,newnode(4)];
enodes(newelem(7),:) = [nbrnode(6),newnode(3),newnode(2)];
enodes(newelem(8),:) = [nbrnode(7),newnode(4),newnode(3)];

M(1,:) = ([Ax,Ay,Az]);
M(2,:) = ([Bx,By,Bz]);
M(3,:) = ([Cx,Cy,Cz]);
M(4,:) = ([Dx,Dy,Dz]);
M(5,:) = ([Ex,Ey,Ez]);

%Springs the region.
i = 1;
imax = 20;
stoptol = 10e-8 * inrad;
tol = stoptol + 10;
while i <= imax & tol > stoptol
    M2(2,1) = sum([(x([(nbrnode([1:3])),badnode]))',M(3,1),M(4,1)]);
    M2(2,2) = sum([(y([(nbrnode([1:3])),badnode]))',M(3,2),M(4,2)]);
    M2(2,3) = sum([(z([(nbrnode([1:3])),badnode]))',M(3,3),M(4,3)]);

    M2(1,1) = sum([(x([(nbrnode([10,11,1]))]))',M(2,1),M(4,1),M(5,1)]);
    M2(1,2) = sum([(y([(nbrnode([10,11,1]))]))',M(2,2),M(4,2),M(5,2)]);
    M2(1,3) = sum([(z([(nbrnode([10,11,1]))]))',M(2,3),M(4,3),M(5,3)]);

    M2(3,1) = sum([(x([(nbrnode([3:6]))]))',M(2,1),M(4,1)]);
    M2(3,2) = sum([(y([(nbrnode([3:6]))]))',M(2,2),M(4,2)]);
    M2(3,3) = sum([(z([(nbrnode([3:6]))]))',M(2,3),M(4,3)]);

    M2(5,1) = sum([(x([(nbrnode([7:10]))]))',M(4,1),M(1,1)]);
    M2(5,2) = sum([(y([(nbrnode([7:10]))]))',M(4,2),M(1,2)]);
    M2(5,3) = sum([(z([(nbrnode([7:10]))]))',M(4,3),M(1,3)]);

    M2(4,1) = sum([(x([(nbrnode([6:7]))]))',M(2,1),M(1,1),M(3,1),M(5,1)]);
    M2(4,2) = sum([(y([(nbrnode([6:7]))]))',M(2,2),M(1,2),M(3,2),M(5,2)]);
    M2(4,3) = sum([(z([(nbrnode([6:7]))]))',M(2,3),M(1,3),M(3,3),M(5,3)]);
    M2 = M2./ 6;
    
    tol = max(sqrt((M2(:,1)-M(:,1)).^2 + (M2(:,2)-M(:,2)).^2));
    
    M = M2;
    x([badnode,([newnode(1:4)])]) = ([M(:,1)]);
    y([badnode,([newnode(1:4)])]) = ([M(:,2)]);
    z([badnode,([newnode(1:4)])]) = ([M(:,3)]);
    i = i + 1;
end

%Creates the output structure.
fem = fem_struct;
fem.e = enodes;
fem.x = x;
fem.y = y;
fem.z = z;
fem.nei = nei; 

%Update the connectivity for the bad node and adds the connections for 
%the new nodes.
fem.nei(badnode,1:nc) = ([newnode([1,3,4]),nbrnode([10,11,1]),zeros(1,nc-6)]);
fem.nei(newnode(1),1:nc) = ([badnode,nbrnode([1:3]),newnode([2,3]),zeros(1,nc-6)]);
fem.nei(newnode(2),1:nc) = ([newnode(1),nbrnode([3:6]),newnode(3),zeros(1,nc-6)]);
fem.nei(newnode(3),1:nc) = ([badnode,newnode([1,2]),nbrnode([6:7]),newnode(4),zeros(1,nc-6)]); %TCM 8/27/04
fem.nei(newnode(4),1:nc) = ([nbrnode([7:10]),badnode,newnode(3),zeros(1,nc-6)]);

%Updates the nei for neighbor nodes. These nodes get an additional
%connector.
spb = ([1,3,6,7,10]);
for i = 1:5
    addconn = ([badnode,newnode([1:4]),badnode]);
    tmp = nei(nbrnode(spb(i)),:);
	ij = find(tmp == badnode);
    ik = min((find(tmp == 0)) + 1);
    if isempty(ik)  % TCM 08/13/2007 -- Added this check in case ik was empty.
       ik= size(nei,2)+1;
    end
    fem.nei(nbrnode(spb(i)),1:ik) = ([tmp(1:(ij-1)),addconn(i+1),addconn(i),tmp((ij+1):(ik-1))]);
end

%Updates the nei for neighbor nodes. These nodes get different connector.
spn = ([2,4,5,8,9]); 
for i = 1:5
    addnode = ([newnode(1),newnode(2),newnode(2),newnode(4),newnode(4)]);
    tmp = nei(nbrnode(spn(i)),:);
    ij = find(tmp == badnode);
    fem.nei(nbrnode(spn(i)),1:length(tmp)) = ([tmp(1:(ij-1)),addnode(i),tmp((ij+1):end)]);
end

%Determine the triqual for the final form to see if has very low quality
%elements.
tmp1 = (nbrelem);
tempL3=(fem.x(fem.e(tmp1,1))-fem.x(fem.e(tmp1,2))).^2+(fem.y(fem.e(tmp1,1))-...
   fem.y(fem.e(tmp1,2))).^2;
tempL1=(fem.x(fem.e(tmp1,2))-fem.x(fem.e(tmp1,3))).^2+(fem.y(fem.e(tmp1,2))-...
   fem.y(fem.e(tmp1,3))).^2;
tempL2=(fem.x(fem.e(tmp1,3))-fem.x(fem.e(tmp1,1))).^2+(fem.y(fem.e(tmp1,3))-...
   fem.y(fem.e(tmp1,1))).^2;
tempLs = ([tempL1 + tempL2 + tempL3]);
xnodes = fem.x(fem.e(tmp1,:));
ynodes = fem.y(fem.e(tmp1,:));
temparea = 0.5*(xnodes(:,1).*(ynodes(:,2)-ynodes(:,3))+xnodes(:,2).*...
   (ynodes(:,3)-ynodes(:,1))+xnodes(:,3).*(ynodes(:,1)-ynodes(:,2)));
fem.ar(tmp1) = temparea;
tempq = (4 * sqrt(3) * temparea) ./ tempLs;
clear tempL3 tempL2 tempL1 tempLs xnodes ynodes temparea;

fem1 = fem;
%Identify very low quality elements and attempts to use a line swap to fix.
poor = find(tempq < .1);
nflag = 1;good = 1;
if ~isempty(poor)
   for it = 1:length(poor)
      je = tmp1(poor(it));
      nodes = fem.e(je,:);
      nflag = 0;
      
      %Determine the angles to see if an edge can be fliped
      a2 = (x(enodes(je,3))-x(enodes(je,2))).^2+(y(enodes(je,3))-y(enodes(je,2))).^2;
      b2 = (x(enodes(je,1))-x(enodes(je,3))).^2+(y(enodes(je,1))-y(enodes(je,3))).^2;
      c2 = (x(enodes(je,2))-x(enodes(je,1))).^2+(y(enodes(je,2))-y(enodes(je,1))).^2;
      A = (180/pi)*acos((b2+c2-a2)./(2*sqrt(b2).*sqrt(c2)));
      B = (180/pi)*acos((c2+a2-b2)./(2*sqrt(c2).*sqrt(a2)));
      C = (180/pi)*acos((a2+b2-c2)./(2*sqrt(a2).*sqrt(b2)));
      [test,ind] = max([A,B,C]);
      if test > 160
         %Find the element numbers to flip the edge.
         nogood = nodes(ind);
         swap = setdiff(nodes,nogood);
         [temp1,tc] = find(enodes == swap(1) | enodes == swap(2));
         [b,j,k] = unique(temp1);
         temp2 = setdiff(1:length(temp1),j);
         swap = temp1(temp2);
         try
            fem1 = line_swap(fem,swap(1),swap(2));
         catch
            fem1 = fem;
         end
         
         %Spring the new mesh after the line swap.
         it = 1;
         while it < 3
         temp = find(fem1.nei(newnode(1),:) ~= 0);
         tempnei = fem1.nei(newnode(1),temp);
         fem1.x(newnode(1)) = mean(fem1.x(tempnei));
         fem1.y(newnode(1)) = mean(fem1.y(tempnei));
         fem1.z(newnode(1)) = mean(fem1.z(tempnei));
         temp = find(fem1.nei(newnode(2),:) ~= 0);
         tempnei = fem1.nei(newnode(2),temp);
         fem1.x(newnode(2)) = mean(fem1.x(tempnei));
         fem1.y(newnode(2)) = mean(fem1.y(tempnei));
         fem1.z(newnode(2)) = mean(fem1.z(tempnei));
         temp = find(fem1.nei(newnode(3),:) ~= 0);
         tempnei = fem1.nei(newnode(3),temp);
         fem1.x(newnode(3)) = mean(fem1.x(tempnei));
         fem1.y(newnode(3)) = mean(fem1.y(tempnei));
         fem1.z(newnode(3)) = mean(fem1.z(tempnei));
         temp = find(fem1.nei(newnode(4),:) ~= 0);
         tempnei = fem1.nei(newnode(4),temp);
         fem1.x(newnode(4)) = mean(fem1.x(tempnei));
         fem1.y(newnode(4)) = mean(fem1.y(tempnei));
         fem1.z(newnode(4)) = mean(fem1.z(tempnei));
         temp = find(fem1.nei(badnode,:) ~= 0);
         tempnei = fem1.nei(badnode,temp);
         fem1.x(badnode) = mean(fem1.x(tempnei));
         fem1.y(badnode) = mean(fem1.y(tempnei));
         fem1.z(badnode) = mean(fem1.z(tempnei));
         it = it + 1;
         end
         
         %Use triqual to determine if the new mesh is better quality.
         tmp1 = nbrelem;
         tempL3=(fem1.x(fem1.e(tmp1,1))-fem1.x(fem1.e(tmp1,2))).^2+(fem1.y(fem1.e(tmp1,1))-...
            fem1.y(fem1.e(tmp1,2))).^2;
         tempL1=(fem1.x(fem1.e(tmp1,2))-fem1.x(fem1.e(tmp1,3))).^2+(fem1.y(fem1.e(tmp1,2))-...
            fem1.y(fem1.e(tmp1,3))).^2;
         tempL2=(fem1.x(fem1.e(tmp1,3))-fem1.x(fem1.e(tmp1,1))).^2+(fem1.y(fem1.e(tmp1,3))-...
            fem1.y(fem1.e(tmp1,1))).^2;
         tempLs = ([tempL1 + tempL2 + tempL3]);
         xnodes = fem1.x(fem1.e(tmp1,:));
         ynodes = fem1.y(fem1.e(tmp1,:));
         temparea = 0.5*(xnodes(:,1).*(ynodes(:,2)-ynodes(:,3))+xnodes(:,2).*...
            (ynodes(:,3)-ynodes(:,1))+xnodes(:,3).*(ynodes(:,1)-ynodes(:,2)));
         fem1.ar(tmp1) = temparea;
         tempq = (4 * sqrt(3) * temparea) ./ tempLs;
         clear tempL3 tempL2 tempL1 tempLs xnodes ynodes temparea;
         if min(tempq) > .4
            fem = fem1;
            nflag = 1;
         else
            nflag = 0;
            good = 6;
         end
      end
   end
end

% If problems with new mesh, then retriangulate

%Correct output if line swap failed to run
if sum(nflag) == 0
   fem = patch_update(fem_struct,fem1,jj);
   good = 6;
end

% Correct output if invalid elements were created.
bad = find(fem.ar < 0);
if ~isempty(bad)
    fem = patch_update(fem_struct,fem1,jj);
    good = 6;
end

%Display message if invalid elements where created.
if good == 6
   disp('The nodally updated mesh contained invalid or badly conditioned');
   disp('elements, therefore the patch was retriangulated which should');
   disp('reduce the connecitivity but is not guaranteed to do so.');
   disp(' ');
end

return 
