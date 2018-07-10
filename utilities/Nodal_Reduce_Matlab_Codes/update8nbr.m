function fem = update8nbr(fem_struct,jj)
% update8nbr takes a node with 8 elements connected to it
% (which defines a region) addes one node and two elements
% to that region then regroups so that no node has more than
% six elements connected inside the region.  The new nodes
% and elements are also sprung.
%
% NOTE:  fem_struct.nei must be a component
%
% CALL:  fem = update8nbr(fem_struct,jj)
%
% Variables
% fem_struct -- is the finite element structure from opnml
% jj -- is the node number that has 8 elements connected to it.
% fem -- is the updated finite element structure
%        Note:  Only fem.x,fem.y,fem.z,fem.e, and fem.nei
%               get updated.  fem.bnd does not need to be updated.
%
%
% Name:   update8nbr.m
% Written By:  Chris Massey
% Date:   July 15, 2003
% Last Modifed:  
%         Sept. 19, 2006 -- Chris Massey, NRL Code 7322, S.S.C. MS 39529
%                      % Add test to make sure the correct nodal neighbor
%                        connectivity existed for the requested update
%                        and added a test to create the nodal neighbor list
%                        if not present.
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
if tempnc ~= 8
    error(['There are ',num2str(tempnc),' nodal neighbors for node number ',...
           num2str(jj),'. This routine only works for 8 nodal neighbors.']);
end

vod = [1:8,1:8];
x = fem_struct.x;
y = fem_struct.y;
z = fem_struct.z;
enodes = fem_struct.e;


[nnodes,nc] = size(nei);
nelems = size(enodes,1);
newndn = nnodes+1; % This is the new node number being added
neweln = nelems+[1,2]; %These are the new element numbers being added

badndn = jj; % This is the bad node number
extndn = nei(jj,1:8); % These are the neighbor nodes of badndn listed in a ccw manner

[nbrelems,J] = find(enodes == badndn); %These are the elements connected to the bad node
clear J

% Finds which set of nodes would be best to add more connectors to.
for i = 1:8
    tem(i) = min([find(nei(extndn(i),:)==0)-1,nc]);
end

temp1 = [max(tem([1,5])),max(tem([2,6])),max(tem([3,7])),max(tem([4,8]))];
temp2 = [min(tem([1,5])),min(tem([2,6])),min(tem([3,7])),min(tem([4,8]))];
k = find(temp1 == min(temp1));
k2 = find(temp2(k) == min(temp2(k)));
j = k(k2(1));

spb = j;
spn = j+4;
clear tem

clear tem temp1 temp2 k k2 j

% First guess of the new coordinates and bathymetry.
Mx = [sum(y([extndn(vod(spb:spn)),badndn])),sum(y([extndn(vod(spn:spn+4)),badndn]))];
My = [sum(x([extndn(vod(spb:spn)),badndn])),sum(x([extndn(vod(spn:spn+4)),badndn]))];
Mz = [sum(z([extndn(vod(spb:spn)),badndn])),sum(z([extndn(vod(spn:spn+4)),badndn]))];
xbar = My/6;
ybar = Mx/6;
zbar = Mz/6;
x([badndn,newndn]) = xbar;
y([badndn,newndn]) = ybar;
z([badndn,newndn]) = zbar;

% Update the effected Elements that need to be updated
[updelem,J] = find(enodes(nbrelems,:) == extndn(vod(spb+5)) | enodes(nbrelems,:) == extndn(vod(spb+7)));
clear J
[updelem2,J] = find(enodes(nbrelems(updelem),:) == badndn);
k = find(enodes(nbrelems(updelem),:) == badndn);
temp = enodes(nbrelems(updelem),:);
temp(k) = newndn;
enodes(nbrelems(updelem),:) = temp;

clear temp J k

% Adds the Two New Elements
enodes(neweln(1),:) = [badndn,newndn,extndn(spb)];
enodes(neweln(2),:) = [extndn(spn),newndn,badndn];

% Springs the Region
imax = 20;
stoptol = 10e-7;
i = 1;
tol = stoptol+10;
while i <= imax & tol > stoptol
    Mx = [sum(y([extndn(vod(spb:spn)),newndn])),sum(y([extndn(vod(spn:spn+4)),badndn]))];
    My = [sum(x([extndn(vod(spb:spn)),newndn])),sum(x([extndn(vod(spn:spn+4)),badndn]))];
    Mz = [sum(z([extndn(vod(spb:spn)),newndn])),sum(z([extndn(vod(spn:spn+4)),badndn]))];
    xbar = My/6;
    ybar = Mx/6;
    zbar = Mz/6;
    tol = norm([x([badndn,newndn])-xbar',y([badndn,newndn])-ybar']);
    x([badndn,newndn]) = xbar;
    y([badndn,newndn]) = ybar;
    z([badndn,newndn]) = zbar;
    i = i + 1;
end

% Creates the Output Structure
fem = fem_struct;
fem.e = enodes;
fem.x = x;
fem.y = y;
fem.z = z;

% Update the Neighbor List for each node
% These correspond to the badndn and the newndn
fem.nei(badndn,1:nc) = [newndn,extndn(vod(spb:spn)),zeros(1,nc-6)];
fem.nei(newndn,1:nc) = [badndn,extndn(vod(spn:spn+4)),zeros(1,nc-6)];

% These are the nodes which get a new connector added two them
temp = nei(extndn(spb),:);
k = min([min(find(temp ==0)),size(temp(:),1)+1]);
j = find(temp == badndn);
fem.nei(extndn(spb),1:k) = [temp(1:j),newndn,temp(j+1:k-1)];
clear temp j k

temp = nei(extndn(spn),:);
k = min([min(find(temp == 0)),size(temp(:),1)+1]);
j = find(temp == badndn);
fem.nei(extndn(spn),1:k) = [temp(1:j-1),newndn,temp(j:k-1)];
clear temp j k

% These are all the other nodes
temp = nei(extndn(vod([spn+1,spn+2,spn+3])),:);
j = find(temp == badndn);
temp(j) = newndn;
fem.nei(extndn(vod([spn+1,spn+2,spn+3])),:) = temp;
clear temp j

%Determine the triqual for the final form to see if has very low quality
%elements.
tmp1 = ([nbrelems;neweln']);
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
poor = find(tempq < .1);

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
         temp = find(fem1.nei(newndn,:) ~= 0);
         tempnei = fem1.nei(newndn,temp);
         fem1.x(newndn) = mean(fem1.x(tempnei));
         fem1.y(newndn) = mean(fem1.y(tempnei));
         fem1.z(newndn) = mean(fem1.z(tempnei));
         temp = find(fem1.nei(badndn,:) ~= 0);
         tempnei = fem1.nei(badndn,temp);
         fem1.x(badndn) = mean(fem1.x(tempnei));
         fem1.y(badndn) = mean(fem1.y(tempnei));
         fem1.z(badndn) = mean(fem1.z(tempnei));

         %Use triqual to determine if the new mesh is better quality.
         tmp1 = ([nbrelems;neweln']);
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
   disp('reduce the connectivity but is not guaranteed to do so.');
   disp(' ');
end


return
