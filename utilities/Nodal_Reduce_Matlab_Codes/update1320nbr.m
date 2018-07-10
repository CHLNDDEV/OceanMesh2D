function fem = update1320nbr(fem_struct,jj)
% update1320nbr takes a node with connectivity between 13 and 20
%  and adds six new nodes and uses delanunay trianglulation to 
%  create the new elements. The region is then sprung. An nei is
%  not required but if one is not present, it will be generated.
%
%  Variables
%  fem_struct -- finite element structure from opnml.
%  jj -- the node number that has 12 elements connected to it.
%  fem -- the updated finite element structure.
%        Note:  Only fem.x,fem.y,fem.z,fem.e, and fem.nei
%               get updated.  fem.bnd does not need to be updated.
%
%  Usage -- fem = update1320nbr(fem_struct,jj)
%
%  FileName: update1320nbr.m
%  Written by: Ben Holladay (SEAP Student 2004)
%  Date: July 7,2004
%  Modified:  August 30, 2004 -- Chris Massey
%             June 14,2005 -- Ben Holladay
%                   Fix: Ensure that delaunay doesn't create elements
%                   outside the patch.
% Last Modifed:  
%         Sept. 19, 2006 -- Chris Massey, NRL Code 7322, S.S.C. MS 39529
%                      % Add test to make sure the correct nodal neighbor
%                        connectivity existed for the requested update
%         July 14, 2008 -- Chris Massey, NRL Code 7322
%               moved test for 4 neighbors to after the test for nei.
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
if (tempnc < 13 | tempnc > 20)
    error(['There are ',num2str(tempnc),' nodal neighbors for node number ',...
           num2str(jj),'. This routine only works for 13-20 nodal neighbors.']);
end

%Sets the intial fem_struct variables.
x = fem_struct.x;
y = fem_struct.y;
z = fem_struct.z;
enodes = fem_struct.e;

%Adds new nodes,determines connectivity, and labels the bad node, 
%neighboring nodesm and neighboring elements.
badnode = jj;
[nbrelem,J] = find(enodes == badnode);
[nnodes,trnc] = size(nei);
nelems = size(enodes,1);
temp1 = sum(ismember(nei(jj,:),0));
nc = trnc - temp1;
nbrnode = nei(badnode,1:nc);
newnode = nnodes + [1:6];
vod = ([1:nc,1:nc]);
inrad = max(sqrt((x(nbrnode)-x(badnode)).^2 + (y(nbrnode)-y(badnode)).^2));

%Computes the center angle of each element.
b2 = (x(nbrnode(:))-x(badnode)).^2 + (y(nbrnode(:))-y(badnode)).^2;
a2 = (x(nbrnode(:))-x(nbrnode([2:nc,1]))).^2 + (y(nbrnode(:))-y(nbrnode([2:nc,1]))).^2;
A = (180/pi)*acos((b2+b2([2:nc,1])-a2)./(2*sqrt(b2).*sqrt(b2([2:nc,1]))));
ang = A;

%Sets the patterns for how to add the nodes.
addnode(1,:) = ([2,2,2,2,2,3]);
addnode(2,:) = ([2,2,3,2,2,3]);
addnode(3,:) = ([2,3,2,3,2,3]);
addnode(4,:) = ([3,3,2,3,3,2]);
addnode(5,:) = ([2,3,3,3,3,3]);
addnode(6,:) = ([3,3,3,3,3,3]);
addnode(7,:) = ([3,3,3,3,3,4]);
addnode(8,:) = ([3,3,4,3,3,4]);

%Determines the appropriate pattern for adding nodes
addnodepat = addnode((nc - 12),:);

%Labels reference nodes to be use to determine the new nodes locations.
spb(1) = 1;
spb(2) = spb(1) + addnodepat(1);
spb(3) = spb(2) + addnodepat(2);
spb(4) = spb(3) + addnodepat(3);
spb(5) = spb(4) + addnodepat(4);
spb(6) = spb(5) + addnodepat(5);

%Determines the total angle that each node will cover.
testang(1) = sum([ang([spb(1):(spb(2)-1)])]);
testang(2) = sum([ang([spb(2):(spb(3)-1)])]);
testang(3) = sum([ang([spb(3):(spb(4)-1)])]);
testang(4) = sum([ang([spb(4):(spb(5)-1)])]);
testang(5) = sum([ang([spb(5):(spb(6)-1)])]);
testang(6) = sum([ang([spb(6):nc])]);

%Test to make sure no node is trying to cover more than half of the 
%circle.
testr = 1;
while testang > 180 & testr ~=0
    hightmp = find(testang > 180);
    spb(vod(hightemp + 1)) = spb(vod(hightemp + 1)) - 1;
    test = find(testang > 180);
	[testr,testc] = size(test);
end    

%First guess of the new coordinates and bathymetry.
M(1,1) = sum(x([(nbrnode([spb(1):spb(2)])),badnode,badnode]));
M(1,2) = sum(y([(nbrnode([spb(1):spb(2)])),badnode,badnode]));
M(1,3) = sum(z([(nbrnode([spb(1):spb(2)])),badnode,badnode]));

M(2,1) = sum(x([(nbrnode([spb(2):spb(3)])),badnode,badnode]));
M(2,2) = sum(y([(nbrnode([spb(2):spb(3)])),badnode,badnode]));
M(2,3) = sum(z([(nbrnode([spb(2):spb(3)])),badnode,badnode]));

M(3,1) = sum(x([(nbrnode([spb(3):spb(4)])),badnode,badnode]));
M(3,2) = sum(y([(nbrnode([spb(3):spb(4)])),badnode,badnode]));
M(3,3) = sum(z([(nbrnode([spb(3):spb(4)])),badnode,badnode]));

M(4,1) = sum(x([(nbrnode([spb(4):spb(5)])),badnode,badnode]));
M(4,2) = sum(y([(nbrnode([spb(4):spb(5)])),badnode,badnode]));
M(4,3) = sum(z([(nbrnode([spb(4):spb(5)])),badnode,badnode]));

M(5,1) = sum(x([(nbrnode([spb(5):spb(6)])),badnode,badnode]));
M(5,2) = sum(y([(nbrnode([spb(5):spb(6)])),badnode,badnode]));
M(5,3) = sum(z([(nbrnode([spb(5):spb(6)])),badnode,badnode]));

M(6,1) = sum(x([(nbrnode([spb(6):nc,1])),badnode,badnode]));
M(6,2) = sum(y([(nbrnode([spb(6):nc,1])),badnode,badnode]));
M(6,3) = sum(z([(nbrnode([spb(6):nc,1])),badnode,badnode]));

M(7,1) = x(badnode);
M(7,2) = y(badnode);
M(7,3) = z(badnode);

M(1,:) = M(1,:) ./ (addnodepat(1) + 3);
M(2,:) = M(2,:) ./ (addnodepat(2) + 3);
M(3,:) = M(3,:) ./ (addnodepat(3) + 3);
M(4,:) = M(4,:) ./ (addnodepat(4) + 3);
M(5,:) = M(5,:) ./ (addnodepat(5) + 3);
M(6,:) = M(6,:) ./ (addnodepat(6) + 3);

x([([newnode(1:6)]),badnode]) = ([M(:,1)]);
y([([newnode(1:6)]),badnode]) = ([M(:,2)]);
z([([newnode(1:6)]),badnode]) = ([M(:,3)]);

%Uses delaunay trianglulation to created a new and updated element list.
tempx = x([badnode,newnode(:)',nbrnode(:)']);
tempy = y([badnode,newnode(:)',nbrnode(:)']);
tri = delaunay(tempx,tempy);
tri2 = tri;

%Convertst the local element list from tri back to its global numberiing.
for i = 1:nc+7
    temp = find(i == tri);
    tmp = ([badnode,newnode(:)',nbrnode(:)']);
    tri2(temp) = tmp(i);
end

%Checks if the elements are in CCW order.
xnodes = x(tri2);
ynodes = y(tri2);
temparea = 0.5*(xnodes(:,1).*(ynodes(:,2)-ynodes(:,3))+xnodes(:,2).*...
    (ynodes(:,3)-ynodes(:,1))+xnodes(:,3).*(ynodes(:,1)-ynodes(:,2)));
Itri = find(temparea < 0);
tri3 = tri2;
tri3(Itri,2) = tri2(Itri,3);
tri3(Itri,3) = tri2(Itri,2);

%Removes elements that would be created outside the patch. 
%!!MODIFICATION!!
oldbnd = detbndy(enodes(nbrelem,:));
newbnd = detbndy(tri3);
oldbnd = sort(oldbnd,2);
newbnd = sort(newbnd,2);
temp11 = find(ismember(newbnd,oldbnd,'rows') ~= 1);
while isempty(temp11) == 0
    newbnd = detbndy(tri3);
    newbnd = sort(newbnd,2);
    temp11 = find(ismember(newbnd,oldbnd,'rows') ~= 1);
    temp12 = newbnd(temp11,:);
    for it = 1:size(temp12,1)
        [tr,tc] = find(tri3 == temp12(it,1) | tri3 == temp12(it,2));
        [ta,tb,tc] = unique(tr);
        temp13 = setdiff(1:length(tr),tb);
        temp14 = tr(temp13);
        tri3(temp14,:) = [];
    end
end    

%Finds and throws out bad elements and adds new elements to the global list. 
tempa = setdiff(1:1:length(enodes),nbrelem);
enodes = enodes(tempa,1:3);
nelems = size(enodes,1);
tri3r = size(tri3,1);
enodes([(nelems + 1):(nelems + tri3r)],:) = tri3;

%Springs the region.
stoptol = 10e-8 * inrad;
tol = stoptol + 10;
i = 1;
imax = 20;
while i <= imax & tol > stoptol
    M2(7,1) = sum([M(2,1),M(1,1),M(3,1),M(5,1),M(4,1),M(6,1)]);
    M2(7,2) = sum([M(2,2),M(1,2),M(3,2),M(5,2),M(4,2),M(6,2)]);
    M2(7,3) = sum([M(2,3),M(1,3),M(3,3),M(5,3),M(4,3),M(6,3)]);

    [r1,c1] = find(enodes == newnode(1));
    nbr1 = unique(enodes(r1,:));
    temp1 = find(nbr1 == newnode(1));
    nbr1(temp1) = [];
	M2(1,1) = sum(x([nbr1]));
    M2(1,2) = sum(y([nbr1]));
    M2(1,3) = sum(z([nbr1]));
    
    [r2,c] = find(enodes == newnode(2));
    nbr2 = unique(enodes(r2,:));
    temp2 = find(nbr2 == newnode(2));
    nbr2(temp2) = [];
	M2(2,1) = sum(x([nbr2]));
    M2(2,2) = sum(y([nbr2]));
    M2(2,3) = sum(z([nbr2]));
    
    [r3,c3] = find(enodes == newnode(3));
    nbr3 = unique(enodes(r3,:));
    temp3 = find(nbr3 == newnode(3));
    nbr3(temp3) = [];
	M2(3,1) = sum(x([nbr3]));
    M2(3,2) = sum(y([nbr3]));
    M2(3,3) = sum(z([nbr3]));
    
    [r4,c4] = find(enodes == newnode(4));
    nbr4 = unique(enodes(r4,:));
    temp4 = find(nbr4 == newnode(4));
    nbr4(temp4) = [];
	M2(4,1) = sum(x([nbr4]));
    M2(4,2) = sum(y([nbr4]));
    M2(4,3) = sum(z([nbr4]));
    
    [r5,c5] = find(enodes == newnode(5));
    nbr5 = unique(enodes(r5,:));
    temp5 = find(nbr5 == newnode(5));
    nbr5(temp5) = [];
	M2(5,1) = sum(x([nbr5]));
    M2(5,2) = sum(y([nbr5]));
    M2(5,3) = sum(z([nbr5]));
    
	[r6,c6] = find(enodes == newnode(6));
    nbr6 = unique(enodes(r6,:));
    temp6 = find(nbr6 == newnode(6));
    nbr6(temp6) = [];
	M2(6,1) = sum(x([nbr6]));
    M2(6,2) = sum(y([nbr6]));
    M2(6,3) = sum(z([nbr6]));
    
    M2(1,:) = M2(1,:) ./ (size((nbr1),1));
    M2(2,:) = M2(2,:) ./ (size((nbr2),1));
    M2(3,:) = M2(3,:) ./ (size((nbr3),1));
    M2(4,:) = M2(4,:) ./ (size((nbr4),1));
    M2(5,:) = M2(5,:) ./ (size((nbr5),1));
    M2(6,:) = M2(6,:) ./ (size((nbr6),1));
    M2(7,:) = M2(7,:) ./ 6;
    
	tol = max(sqrt((M2(:,1)-M(:,1)).^2 + (M2(:,2)-M(:,2)).^2));
    
    M = M2;;
    x([([newnode(1:6)]),badnode]) = ([M(:,1)]);
    y([([newnode(1:6)]),badnode]) = ([M(:,2)]);
    z([([newnode(1:6)]),badnode]) = ([M(:,3)]);
    i = i + 1;
end    

%Updates the nei for the region.
tempnode = ([newnode,badnode,nbrnode]);
tempcount = size(tempnode,2);
for i = 1:tempcount 
    [tr,tc] = find(tempnode(i) == enodes);
    tempnbr = unique(enodes(tr,:));
    tmp1 = find(tempnbr == tempnode(i));
    tempnbr(tmp1) = [];
    tempnbrx = x(tempnbr) - x(tempnode(i));
    tempnbry = y(tempnbr) - y(tempnode(i));
    angle = (atan2(tempnbry,tempnbrx))*(180/pi);
    tempang = find(angle < 0);
    angle(tempang) = angle(tempang) + 360;
    [tmp2,tempi] = sort(angle);
    tempnbr1 = tempnbr(tempi);
    tmpr = size(tempnbr,1);
    nei(tempnode(i),:) = ([tempnbr1(1:tmpr)',zeros(1,trnc-tmpr)]);     
end    

%Createst the output structure.
fem = fem_struct;
fem.x = x;
fem.y = y;
fem.z = z;
fem.e = enodes;
fem.nei = nei;
return
