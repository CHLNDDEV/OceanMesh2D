function fem = line_swap(fem_struct,j,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Line_swap is mesh refining tool that swaps an edge that two elements
% share. It reconnects the elements so that the edge they share has been 
% flipped. This routine only operates on one line at a time and can only be
% used if the two elements share an edge. This routine can operate on any
% elemnet, even boundry elements, unless the new elements formed by the line
% swap would overlap existing elements. No nodes are moved by this routine.
% The fem.e and fem.ar are the only fields that are updated. If no .nei file
% is present then none will be generated; however, if one is present then
% it will be updated to reflect the new connectivity.
%
% Calls: is_valid_struct.m
%
% Usage: fem = line_swap(fem_struct,j,k);
%
% Variables:
%  fem -- the new split finite element mesh.
%  fem_struct -- the finite element grid structure from the opnml suite.
%          Note: fem_struct.nei will not be generated if not present
%  j -- element number for one element used in the line swap.
%  k -- element number for one element used in the line swap.
%          Note: j and k can be interchanged.
%
% Filename: line_swap.m
% Created by: Ben Holladay
% Date: June 8, 2005
% Last Modified:
%       Oct. 26, 2007 -- Chris Massey, NRL 7322
%       Changed when line swap can not be performed due to mesh tangling,
%       from an error call to a display message and return original fem.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ensure the apropriate number of input arguements is given, checks that the 
%fem_struct is valid.
if nargin == 0 
    error('Not enough input arguments; need a fem_struct and two element numbers.');
end
if ~is_valid_struct(fem_struct)
   error('Input argument to line swap must be a valid fem_struct.');
end
if nargin == 1
    error(['Input arguments must include a valid fem_struct and two element ' , ...
    'numbers that share an edge.']);
elseif nargin == 2
    error('Two element numbers must be given.');
elseif nargin == 3
    if j <= 0 | j > size(fem_struct.e,1) | j ~= floor(j)
        error(['Element numbers must positive integers that are less than ',... 
            'the maximum number of elements.']);
    end
    if k <= 0 | k > size(fem_struct.e,1) | k ~= floor(k)
        error(['Element numbers must positive integers that are less than ',... 
            'the maximum number of elements.']);
    end
else
    error('Too many input arguments.');
end


%Sets the intial fem_struct variables.
x = fem_struct.x;
y = fem_struct.y;
enodes = fem_struct.e;
ar = fem_struct.ar;

%Determine if the elements share a common edge.
elems = enodes([j;k],:);
comp = unique(intersect(elems(1,:),elems(2,:)));
if length(comp) < 2
    error('The elements do not share a common edge.');
end

%Indentify each node for swapping.
[nodes,ti,tj] = unique(elems);
temp1 = setdiff((1:6),ti);
temp1 = elems(temp1);
temp2 = setdiff(nodes,temp1);

%Determine if the angle between the elements is too big for line_swap to
%effectively operate on it.
a2 = (x(enodes([j;k],3))-x(enodes([j;k],2))).^2+(y(enodes([j;k],3))-y(enodes([j;k],2))).^2;
b2 = (x(enodes([j;k],1))-x(enodes([j;k],3))).^2+(y(enodes([j;k],1))-y(enodes([j;k],3))).^2;
c2 = (x(enodes([j;k],2))-x(enodes([j;k],1))).^2+(y(enodes([j;k],2))-y(enodes([j;k],1))).^2;
A = (180/pi)*acos((b2+c2-a2)./(2*sqrt(b2).*sqrt(c2)));
B = (180/pi)*acos((c2+a2-b2)./(2*sqrt(c2).*sqrt(a2)));
C = (180/pi)*acos((a2+b2-c2)./(2*sqrt(a2).*sqrt(b2)));
ang = [A,B,C];
temp3 = find(temp1(1) == elems);
temp4 = find(temp1(2) == elems);
temp5 = ang(temp3);
temp6 = ang(temp4);
temp5 = sum(temp5);
temp6 = sum(temp6);
if temp5 >= 180 | temp6 >= 180
    disp(['The line cannot be swapped because the new elements would '...
        'create overlapping elements.']);
    fem = fem_struct;
    return
elseif temp5 >= 135 | temp6 >= 135
    disp(['Warning: The elements created will contain large angles.']);
end
clear a2 b2 c2 A B C ang

%Swap the lines in the enodes list.
temp3 = setdiff(elems(1,:),temp1);
temp4 = setdiff(elems(2,:),temp1);
temp5 = find(elems(1,:) == temp1(1));
elems(1,temp5) = temp4;
temp6 = find(elems(2,:) == temp1(2));
elems(2,temp6) = temp3;
clear temp5 temp6;

%Recomputes the areas for the new elements.
xnodes = x(elems);
ynodes = y(elems);
temparea = 0.5*(xnodes(:,1).*(ynodes(:,2)-ynodes(:,3))+xnodes(:,2).*...
(ynodes(:,3)-ynodes(:,1))+xnodes(:,3).*(ynodes(:,1)-ynodes(:,2)));
clear xnodes ynodes x y;

%Addes the new elements and areas to the global lists.
enodes([j;k,],:) = elems;
ar([j;k]) = temparea(:);
clear elems temparea;

%Create the output structure.
fem = fem_struct;
fem.e = enodes;
fem.ar = ar;

%Updates the nei if present.
try
    nei = fem_struct.nei;
    nei = [nei,zeros((size(nei,1)),1)];
    tempnei = nei([temp1(:);temp2(:)],:);
    temp3 = find(tempnei(1,:) == temp1(2));
    temp4 = setdiff(1:(size(tempnei,2)),temp3);
    tempnei(1,:) = ([tempnei(1,temp4),0]);
    temp3 = find(tempnei(2,:) == temp1(1));
    temp4 = setdiff(1:(size(tempnei,2)),temp3);
    tempnei(2,:) = ([tempnei(2,temp4),0]);
    temp3 = find(tempnei(3,:) == temp1(1) | tempnei(3,:) == temp1(2));
    if abs(temp3(1)-temp3(2)) == 1
		tempnei(3,:) = [tempnei(3,(1:temp3(1))),temp2(2),tempnei(3,((temp3(2)):(end-1)))];
    else
        tempnei(3,:) = [temp2(2),tempnei(3,(1:end-1))];
    end
    temp3 = find(tempnei(4,:) == temp1(1) | tempnei(4,:) == temp1(2));
    if abs(temp3(1)-temp3(2)) == 1    
        tempnei(4,:) = [tempnei(4,(1:temp3(1))),temp2(1),tempnei(4,((temp3(2)):(end-1)))];
    else
        tempnei(4,:) = [temp2(1),tempnei(4,(1:end-1))];
    end
    nei([temp1(:);temp2(:)],:) = tempnei;
    temp3 = find(nei(:,end) ~= 0);    
    if isempty(temp3) == 1
        nei = nei(:,(1:end-1));
    end
    temp3 = find(nei(:,end) ~= 0);
    if isempty(temp3) == 1
        nei = nei(:,(1:end-1));
    end
    fem.nei = nei;
catch
end

return