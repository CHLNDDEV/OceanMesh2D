function nei = ele2nei(enodes,x,y);
% ele2nei - build a FEM neighbor list from element and node lists
%  ele2nei computes the neighbor list for a FEM mesh specified
%  by the triangular element lists.  Each node's neighbor list
%  is then sorted into counter-clockwise order.
%
%  The resulting neighbor list can be passed to WRITE_NEI
%  to output a FEM neighbor file (.nei) to disk.
%
%  This is a rewrite of the origial ele2nei.m written by
%  Brian O. Blanton in March 1996 which used a .mex file
%  written in C.  It is much faster than the original version
%  and no longer has a maximum of 20 nodal neighbors.
%     
% INPUTS:   e - 3-column element connectivity list (REQ)
%           x - x-coordinate list (REQ)
%           y - y-coordinate list (REQ)
% OUTPUTS:  nei - neighbor list (REQ)
%
%    Use: nei = ele2nei(e,x,y);
%
%
%   Written By:  Ben Holladay
%   Date:        July, 2005
%   Last Modified:
%        Aug. 3, 2010 -- Chris Massey, USACE-ERDC-CHL, Vicksburg, MS
%                Cleaned up the code by clearing variables once they were
%                no longer needed.  This allows for a better use of memory
%                which can be a problem with really large problems.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine number of nodes.
nnodes = size(x,1);

%Determine the maximum nodal connectivity.
tempenodes = enodes(:);
szenodes = size(enodes);
[temp,in] = sort(tempenodes);
temp2 = [temp(2:end);temp(1)];
pbend = temp2-temp;
test = find(pbend ~= 0);
test2 = [0;test(1:end-1)];
pblong = max(test-test2);
clear test test2 temp2
clear temp in pbend

%Create a zeros matrix for nei
nei = zeros(nnodes,pblong*2);

%Begin loop to create the nei. The loop uses the unique command to find the
%element a node is a part of. It then addes the other nodes in this element
%to that nodes nei. The loop runs a number of times equal to the maximum
%connectivity. After each iteration it replaces each occurance of a node
%that was found with zeros to ensure it doesn't find the same nodes twice.
for it = 1:pblong
   %Find all unique nodes still in the list.
   [b,j,k] = unique(tempenodes);
   temp3 = find(b == 0);
   
   %Remove zeros from the unique list.
   j(temp3) = [];
   b(temp3) = [];
   
   %Find the element each node is a part of.
   [inr,inc] = ind2sub(szenodes,j);
   
   %Find the other nodes in each element that a node is a part of.
   temp4 = zeros(size(b,1),2);
   temp5 = find(inc == 1);
   temp6 = repmat([2,3],length(temp5),1);
   temp4(temp5,:) = temp6;
   temp5 = find(inc == 2);
   temp6 = repmat([1,3],length(temp5),1);
   temp4(temp5,:) = temp6;
   temp5 = find(inc == 3);
   temp6 = repmat([1,2],length(temp5),1);
   temp4(temp5,:) = temp6;
   
   %Determine the position in the original enodes list of a node's
   %neighbors to be added to the nei.
   temp7 = sub2ind(szenodes,inr,temp4(:,1));
   temp8 = sub2ind(szenodes,inr,temp4(:,2));
   
   %Adds the neighbor nodes to the nei.
   nei(b,[2*it-1]) = enodes(temp7);
   nei(b,[2*it]) = enodes(temp8);
   
   %Sets all the nodes found in this pass to zero so they will not be
   %found in the next pass.
   tempenodes(j) = 0;
   clear temp3 temp4 temp5 temp6 temp7 temp8
   clear inr inc
end

clear tempenodes enodes
clear j b k


%Removes nodes all duplicate nodes in the nei.
snei = sort(-1*nei,2) * -1;
clear nei
snei = [snei,zeros(size(snei,1),1)];
tsnei = snei(:,2:end);
dupnei = snei(:,1:end-1) - tsnei;
clear tsnei
temp9 = find(dupnei == 0);
snei(temp9) = 0;
s2nei = sort(-1*snei,2) * -1;
clear snei temp9

%Removes the extra zeros in the nei.
[tr,tc] = find(s2nei ~= 0);
nc = max(tc);
wnei = s2nei(:,1:nc);
clear s2nei
%Sorts the nodes into CCW using atan2.

%Ensure zeros in nei can be indexed by atan2.
temp = find(wnei == 0);
wnei(temp) = nnodes+1;
x(nnodes+1) = 0;
y(nnodes+1) = 0;

%Determines the angle for sorting.
temp1 = repmat((1:nnodes)',1,nc);
wrong = atan2(y(wnei)-y(temp1),x(wnei)-x(temp1));

%Ensure zeros in the nei are placed at the end.
wrong(temp) = wrong(temp) + 100;
wnei(temp) = 0;
clear temp

%Change negative angles to positive angles.
%Ensure ordering is the same as the original.
temp2 = find(wrong < 0);
wrong(temp2) = wrong(temp2) + 2*pi;

%Sort the angles and use to the index to order the nodes.
[temp2,right] = sort(wrong,2);
true = sub2ind(size(wnei),temp1,right);
nei = wnei(true);

return
