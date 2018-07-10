function fem = reduce_con_int(fem_struct,nvl)
% reduce_con_int updates all interior nodes that are connected
%  to more than 7 nodes and less than 21 so that every node is
%  connected to at most 7 nodes. A springing routine is used locally
%  to help smooth out the updates.  Note that specifying NVL one can 
%  fix nodes with connectivity of NVL or higher and leave the others.
%
%  NOTE:  fem_struct.nei should already be a component, if not one is
%         computed.
%
%
%  Variables
%  fem_struct -- is the finite element structure from the opnml suite.
%  fem -- the updated finite element structure.
%        Note:  fem.x,fem.y,fem.z,fem.e, fem.nei, and fem.ar
%                  get updated.
%               fem.bnd does not need to be updated.
%  nvl -- lower limit of connectivity to fix should be 8 <= nvl <= 20
%
%  Usage -- fem = reduce_con_int(fem_struct);
%
%  Name: reduce_con_int.m
%  Written by: Ben Holladay (SEAP Student 2004)
%  Date: June 22,2004
%  Modified:  Aug. 30, 2004, Chris Massey
%             Jan. 10, 2006, Chris Massey
%             Aug. 15, 2006, Chris Massey  - Cleaned up comments.
%             Sept. 19, 2006, Chris Massey -- fixed bug
%             Aug. 22, 2007, Chris Massey -- fixed bug
%             Apr. 25, 2018, Chris Massey -- Added nvl option to fix nodes
%                with connectivity higher than or equal to nvl.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ensure the appropriate number of input arguements is given, checks that the 
%fem_struct is valid.
if nargin == 0 
    error('Not enough input arguments; need a valid fem_struct.');
end
if ~is_valid_struct(fem_struct)
   error('Input argument to reduce_con_bnd must be a valid fem_struct.');
end

%Make sure that the lower limit of connecitivity count falls between
% 8 and 20 nodes
if (exist('nvl','var'))
    if nvl < 8
        disp('Connectivity Lower Limit must be at least 8.  Resetting');
        nvl = 8;
    end
    if nvl > 20
        disp('Connectivity Lower Limit must not be greater than 20. Resetting');
        nvl = 20;
    end
else
    nvl = 8;
end

%Sets the intial fem_struct variables and determines the max
%connectivity. Displays errors if the connectivity is to high.

% List of all nodes not on the boundary;
tmpnodes = setdiff(1:1:length(fem_struct.x),unique(fem_struct.bnd(:)));
bndrynodes = unique(fem_struct.bnd(:));

[nelems,~] = size(fem_struct.e);

% Determine if fem_struct.nei is present, if not added it.
try
    [~,nc] = size(fem_struct.nei(tmpnodes,:));
catch
    disp(' ');
    disp('A neighbor list was not present in the finite element structure.');
    disp('One is being added now.');
    fem_struct.nei = ele2nei(fem_struct.e,fem_struct.x,fem_struct.y);
    disp(' ');
    disp('The neigbor list was successfully added.');
    [~,nc] = size(fem_struct.nei(tmpnodes,:));
end

nelems_orig = nelems;
nnodes_orig = length(fem_struct.x);
nc_orig = nc;

[~,J1] = find(fem_struct.nei(bndrynodes,:)~=0);
highbndry = max(J1);
clear J1


if nc > 21
    error(['The connectivity is too high and must be brought down to at least',...
        '20 by hand.']);
end
if nc <= nvl-1 %7
    disp(['The grids nodal connectivity is not higher than ',num2str(nvl),', nothing will be done.']);
    disp('Returning the original mesh.');
    fem = fem_struct;
    return
end

disp(' ')
disp('NOTE: NO UPDATING IS DONE FOR BOUNDARY NODES.');
disp(' ')

%Begins loop to update nodes.
tempr = size(fem_struct.nei,2);
if tempr < nvl-1 %7
    tempr = 0;
else
    tmpr = length(find(fem_struct.nei(tmpnodes,nvl)~=0));
end
imax = 20;i = 1;
i2flag = 1;
fixedint=[];
disp(['Beginning the update loop'])
disp(' ')
while tempr ~= 0 && i <= imax
    disp(['Pass number ',num2str(i,'%4.0f'),' through the loop'])
    
    %tmpnodes2 = tmpnodes;
    %iupdnn=find(fem_struct.nei(tmpnodes,nvl)~=0);  %nodes with high connectivity
    %tmpnodes = tmpnodes(iupdnn);
    %clear tmpnodes2
    [nnodes,~] = size(fem_struct.nei(tmpnodes,:));
    
    j2 = 1;
    while j2 <= nnodes 
        j = tmpnodes(j2);
        temp2 = fem_struct.nei(j,nvl); %8;
        if  temp2 ~= 0
            testnei = fem_struct.nei(j,:);
            % TCM -- Begin -- 09/19/2006
            %temp1 = sum(ismember((fem_struct.nei(j,:)),0)); TCM
            %tempnc = nc - temp1; TCM
            tempnc = sum(~ismember((fem_struct.nei(j,:)),0)); %TCM
            % TCM -- END -- 09/19/2006
            
            disp(['Updating node number ',num2str(j),' using',' fix ',num2str(tempnc)]);
            switch tempnc
                case {13,14,15,16,17,18,19,20}
                    fem_struct = update1320nbr(fem_struct,j);
                case 12
                    fem_struct = update12nbr(fem_struct,j);
                case 11
                    fem_struct = update11nbr(fem_struct,j);
                case 10
                    fem_struct = update10nbr(fem_struct,j);
                case 9
                    fem_struct = update9nbr(fem_struct,j);
                case 8
                    fem_struct = update8nbr(fem_struct,j);
                otherwise
                    disp('Already fixed.')
            end
            if sum(testnei) == sum(fem_struct.nei(j,:)) %?? What is this test?
                fixedint(i2flag) = j;
                temp1 = find(fem_struct.nei(j,:) ~= 0,1,'last');
                fixedintnei(i2flag) = temp1;
                i2flag = i2flag + 1;
            end
        end
        if size(fem_struct.nei,2) < nvl %8
            j2 = nnodes + 1;
        else
            j2 = j2 + 1;
        end
    end
    
    %Check if any nodes where added to the do not update lists.
    try
        fixedint = fixedint;
        fixedintnei = fixedintnei;
    catch
        fixedint = 1:0;
        fixedintnei = 1:0;
    end
    
    %Create do not update list for interior nodes.
    
    %newnei = fem_struct.nei(fixedint,:);    
    %temp1 = ismember(fem_struct.nei(fixedint,:),0);
    %temp2 = sum(temp1,2);
    %temp3 = [(size(fem_struct.nei,2))-temp2];
    %temp4 = fixedintnei' - temp3;
    %temp5 = find(temp4 < 0);
    %fixedint(temp5) = [];
    %fixedintnei(temp5) = [];
    %temp6 = size(temp5,1);
    %i2flag = i2flag - temp6;

    int = setdiff(1:length(fem_struct.x),bndrynodes);
    
    tmpnodes = setdiff(int,fixedint);
    
    i = i + 1;
    %tmpnodes = setdiff(1:1:length(fem_struct.x),unique(fem_struct.bnd(:)));
    if size(fem_struct.nei,2) > nvl-1
       tempr = length(find(fem_struct.nei(tmpnodes,nvl) ~= 0));
    else
        tempr = 0;
    end
end



%Creates the output structure.
fem = fem_struct;
[~,J1] = find(fem_struct.nei(tmpnodes,:)~=0);
highint = max(J1);

fem.nei = fem.nei(:,1:max(highint,highbndry));
fem = el_areas(fem);


% Display a message that we reached the maximum number of loop iterations
% prior to minimizing the interior nodal connectivity list
if i>imax && highint>nvl-1 %7
    disp(' ')
    disp(['Maximum loop iteration count of ',num2str(imax),' was',...
        ' exceeded before all interior nodes']);
    disp('reached their miniminal connecitvity.');
end

% Display a summary of mesh updates.
disp(' ')
disp([num2str(size(fem.e,1)-nelems_orig),' elements have been added']);
disp([num2str(length(fem.x)-nnodes_orig),' nodes have been added']);
disp(['Maximum interior nodal connectivity was reduced from ',...
    num2str(nc_orig), ' to ',num2str(highint),'.']);
disp(['The maximum nodal connectivity for a boundary node',...
    ' is ',num2str(highbndry),'.']);
disp(' ');

return
