function [fem_grid_struct,bndry_struct]=read_adcirc_mesh(fort14name)
%READ_ADCIRC_MESH converts an ADCIRC grd file to an OPNML fem_grid_struct
% and a boundary_struct.
%
% ADCIRC grid information assumed in "fort.14" format.
%
% Input:  fort14name - path/name of fort.14 file;  if not passed,
%                      assumes fort.14 in the currect working dir.
% Output: fem_grid_struct - OPNML grid structure
%         bndry_struct -- ADCIRC Boundary Information in a structure
%
% Usage:
%       [fem_grid_struct,bndry_struct]=read_adcirc_mesh(fort14name);
%
% Calls: opnml
%
% This is a modified version of the file GRD_TO_OPNML that ignores
% the comment lines in a fort.14 file and reads in the boundary
% information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Modified by:  Chris Massey, January 13, 2003.
%                    Chris Massey, NRL Code 7322 -- Aug. 8, 2007
%                    -- Added boundry structure read/output.
%      Aug. 14, 2007 -- Chris Massey, NRL Code 7322
%            Added a test to catch the case when there is no boundary
%            information in the file and a test to see when to read in
%            just the node and element list based on the number of
%            requested outputs.
%      July 1, 2008 -- Chris Massey
%            Added a test to catch the cases when there are no open
%            or no closed boundary types.
%      July 2, 2008 -- Chris Massey, NRL Code 7322
%            Added conditions to handle when comments are present or not
%            in the boundary file.  Added Weir and Pipe type to the flow
%            boundary types.  Renamed and expanded bndry_struct components.
%      July 2, 2008 -- Chris Massey
%            Accounted for the case when one asks for a boundary_struct
%            and it is not present in file.
%            Accounted for the case when there are only 1 or 2 nodes
%            per section.
%      Aug. 3, 2010 -- Chris Massey, USACE-ERDC-CHL, Vicksburg, MS
%            Expanded to include the field values for weir, weir pair,
%            and pipe pairs in the structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<1
    % assume fort.14 filename in the current wd.
    fort14name='fort.14';
end

if nargout <= 1
    disp('Only the node and element list will be read in from the file.');
end

% Open fort.14 file
[f14,message]=fopen(fort14name,'r');
if (f14<0)
    error(message)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get grid info
gridname = fscanf(f14,'%s',1); % Reads in Grid_Name
fgets(f14); %moves the pointer to the next line (skipping any comments).

% Read in the Number of Elements and Number of Nodes
temp=fscanf(f14,'%d %d',2);
nn=temp(2);
ne=temp(1);
fgets(f14); % moves the pointer to the next line (skipping any comments).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get node locations
temp1 = fscanf(f14,'%d %f %f %f',[4 1])'; % Reads in First set of Nodes
x(1,1) = temp1(2);
y(1,1) = temp1(3);
z(1,1) = temp1(4);
fgets(f14); %moves the pointer to the next line (skipping any comments).

temp=fscanf(f14,'%d %f %f %f',[4 nn-1])'; % Reads in all other Nodes
x(2:nn,1)=temp(:,2);
y(2:nn,1)=temp(:,3);
z(2:nn,1)=temp(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get elements
temp=fscanf(f14,'%d %d %d %d %d',[5 1])'; % Reads in the First Element
e(1,1:3)=temp(:,3:5);
fgets(f14); %moves the pointer to the next line (skipping any comments).

% Read in the Remaining Elements
temp=fscanf(f14,'%d %d %d %d %d',[5 ne-1])';
e(2:ne,1:3)=temp(:,3:5);
if ~feof(f14)
    fgets(f14);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fem_grid_struct.name=gridname;
fem_grid_struct.x=x;
fem_grid_struct.y=y;
fem_grid_struct.z=z;
fem_grid_struct.e=e;
fem_grid_struct.bnd=detbndy(e);
fem_grid_struct=el_areas(fem_grid_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section reads in the boundary information and stores it in a
% a boundary structure, (bndry_struct);

if (nargout == 2 & ~feof(f14))
    % Read number of elevation boundary segments
    pos = ftell(f14);
    data = fgets(f14);
    if isempty(str2num(data))
        fseek(f14,pos,'bof');
        data = fscanf(f14,'%d\n');
        fgets(f14); %Alpha numeric junk
    else
        %No alpha numeric junk present
        data = str2num(data);
    end
    nelevbndy = data(1);%fscanf(f14,'%d\n',1)
    %fgets(f14) %Alpha numeric junk
    nelevbndy;
    if isempty(nelevbndy)
        disp('No boundary information was found in the file.')
        bndry_struct = [];
    else
        pos = ftell(f14);
        data = fgets(f14);
        if isempty(str2num(data))
            fseek(f14,pos,'bof');
            data = fscanf(f14,'%d\n');
            fgets(f14); %Alpha numeric junk
        else
            %No alpha numeric junk present
            data = str2num(data);
        end
        nnodes_elev = data(1);
        pb = 1;
        for i=1:nelevbndy
            pos = ftell(f14);
            data = fgets(f14);
            if isempty(str2num(data))
                fseek(f14,pos,'bof');
                data = fscanf(f14,'%d\n');
                fgets(f14); %Alpha numeric junk
            else
                %No alpha numeric junk present
                data = str2num(data);
            end
            elev_seg(i) = data(1);
            np = data(1);
            if length(data)>1
                elev_type(i) = data(2);
            else
                elev_type(i) = 0;
            end
            %Read the first line
            pos = ftell(f14);
            data = fgets(f14);
            if isempty(str2num(data))
                fseek(f14,pos,'bof');
                data = fscanf(f14,'%d\n',1);
                fgets(f14); %Alpha numeric junk
            else
                %No alpha numeric junk present
                data = str2num(data);
            end
            elev_nodes(pb) = data(1);
            if (np > 2)
                %Read lines 2 through np-1 ; usually no comments
                data = fscanf(f14,'%d\n',[1,np-2]);
                elev_nodes(pb+1:pb+np-2) = data(1:np-2);
            end
            if(np > 1)
                %Read line np
                pos = ftell(f14);
                data = fgets(f14);
                if isempty(str2num(data))
                    fseek(f14,pos,'bof');
                    data = fscanf(f14,'%d\n',1);
                    fgets(f14); %Alpha numeric junk
                else
                    %No alpha numeric junk present
                    data = str2num(data);
                end
                elev_nodes(pb+np-1) = data(1);
            end
            %Advance the counter
            pb = pb + np;
            clear data
        end
        if (nelevbndy==0)
            elev_seg = [];
            elev_nodes = [];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Process flow boundary nodes
        pos = ftell(f14);
        data = fgets(f14);
        if isempty(str2num(data))
            fseek(f14,pos,'bof');
            data = fscanf(f14,'%d\n');
            fgets(f14); %Alpha numeric junk
        else
            %No alpha numeric junk present
            data = str2num(data);
        end
        nflowbndy = data(1);
        %Read the first line
        pos = ftell(f14);
        data = fgets(f14);
        if isempty(str2num(data))
            fseek(f14,pos,'bof');
            data = fscanf(f14,'%d\n');
            fgets(f14); %Alpha numeric junk
        else
            %No alpha numeric junk present
            data = str2num(data);
        end
        nnodes_flow = data(1);

        pbl = 1;
        pbw = 1;
        pbwp = 1;
        pbpp = 1;
        flow_seg = [];
        flow_type = [];
        flow_nodes = [];
        weir_nodes = [];
        weir_pair_nodes = [];
        pipe_pair_nodes = [];
        weir_values = [];
        weir_pair_values = [];
%         weir_barlanht = [];
%         weir_barlancfsp = [];
%         weir_pair_barinht = [];
%         weir_pair_barincfsb = [];
%         weir_pair_barincfsp = [];
        pipe_pair_values = [];

        % for all these need to put in the fgets test as above to
        % remove comments/no comments effects

        for i=1:nflowbndy
            pos = ftell(f14);
            data = fgets(f14);
            if isempty(str2num(data))
                fseek(f14,pos,'bof');
                data = fscanf(f14,'%d %d',[2,1]);
                fgets(f14); %Alpha numeric junk
            else
                %No alpha numeric junk present
                data = str2num(data);
            end
            flow_seg(i) = data(1);
            flow_type(i) = data(2);
            np = data(1);
            clear data
            %fgets(f14); %Alpha numeric junk
            if (ismember(flow_type(i),[3,4,5,13,23,24,25])==0)
                %Read the first line
                pos = ftell(f14);
                data = fgets(f14);
                if isempty(str2num(data))
                    fseek(f14,pos,'bof');
                    data = fscanf(f14,'%d\n',1);
                    fgets(f14); %Alpha numeric junk
                else
                    %No alpha numeric junk present
                    data = str2num(data);
                end
                flow_nodes(pbl) = data(1);
                if(np > 2)
                    %Read lines 2 through np-1 ; usually no comments
                    data = fscanf(f14,'%d\n',[1,np-2]);
                    flow_nodes(pbl+1:pbl+np-2) = data(1:np-2);
                end
                if(np > 1)
                    %Read line np
                    pos = ftell(f14);
                    data = fgets(f14);
                    if isempty(str2num(data))
                        fseek(f14,pos,'bof');
                        data = fscanf(f14,'%d\n',1);
                        fgets(f14); %Alpha numeric junk
                    else
                        %No alpha numeric junk present
                        data = str2num(data);
                    end
                    flow_nodes(pbl+np-1) = data(1);
                end
                %Advance the counter
                pbl = pbl + np;

            elseif (ismember(flow_type(i),[3,13,23])==1)
                %Read the first line
                pos = ftell(f14);
                data = fgets(f14);
                if isempty(str2num(data))
                    fseek(f14,pos,'bof');
                    data = fscanf(f14,'%d %f %f\n',[3,1]);
                    fgets(f14); %Alpha numeric junk
                else
                    %No alpha numeric junk present
                    data = str2num(data);
                end
                weir_nodes(pbw) = data(1);
                weir_values(pbw,1:2) = data(2:3)';
                if(np > 2)
                    %Read lines 2 through np-1 ; usually no comments
                    data = fscanf(f14,'%d %f %f\n',[3,np-2]);
                    weir_nodes(pbw+1:pbw+np-2) = data(1,1:np-2);
                    weir_values(pbw+1:pbw+np-2,1:2) = data(2:3,1:np-2)';
                end
                if(np > 1)
                    %Read line np
                    pos = ftell(f14);
                    data = fgets(f14);
                    if isempty(str2num(data))
                        fseek(f14,pos,'bof');
                        data = fscanf(f14,'%d %f %f\n',[3,1]);
                        fgets(f14); %Alpha numeric junk
                    else
                        %No alpha numeric junk present
                        data = str2num(data);
                    end
                    weir_nodes(pbw+np-1) = data(1);
                    weir_values(pbw+np-1,1:2) = data(2:3)';
                end
                %Advance the counter
                pbw = pbw + np;
            elseif(ismember(flow_type(i),[4,24])==1)
                % Read the first line
                pos = ftell(f14);
                data = fgets(f14);
                if isempty(str2num(data))
                    fseek(f14,pos,'bof');
                    data = fscanf(f14,'%d %d %f %f %f\n',[5,1]);
                    fgets(f14); %Alpha numeric junk
                else
                    %No alpha numeric junk present
                    data = str2num(data);
                end
                weir_pair_nodes(pbwp,1:2) = data(1:2);
                weir_pair_values(pbwp,1:3) = data(3:5);
                
                if(np > 2)
                    %Read lines 2 through np-1 ; usually no comments
                    data = fscanf(f14,'%d %d %f %f %f\n',[5,np-2]);
                    weir_pair_nodes(pbwp+1:pbwp+np-2,1:2) = data(1:2,1:np-2)';
                    weir_pair_values(pbwp+1:pbwp+np-2,1:3) = data(3:5,1:np-2)';
                end
                if(np > 1)
                    %Read line np
                    pos = ftell(f14);
                    data = fgets(f14);
                    if isempty(str2num(data))
                        fseek(f14,pos,'bof');
                        data = fscanf(f14,'%d %d %f %f %f\n',[5,1]);
                        fgets(f14); %Alpha numeric junk
                    else
                        %No alpha numeric junk present
                        data = str2num(data);
                    end
                    weir_pair_nodes(pbwp+np-1,1:2) = data(1:2);
                    weir_pair_values(pbwp+np-1,1:3) = data(3:5)';
                end
                %Advance the counter
                pbwp = pbwp + np;

            elseif(ismember(flow_type(i),[5,25])==1)
                %Read first line for this section
                pos = ftell(f14);
                data = fgets(f14);
                if isempty(str2num(data))
                    fseek(f14,pos,'bof');
                    data = fscanf(f14,'%d %d %f %f %f %f %f %f\n',[8,1]);
                    fgets(f14); %Alpha numeric junk
                else
                    %No alpha numeric junk present
                    data = str2num(data);
                end
                pipe_pair_nodes(pbpp,1:2) = data(1:2);
                pipe_pair_values(pbpp,1:6) = data(3:8);
                if(np > 2)
                    %Read lines 2 through np-1 ; usually no comments
                    data = fscanf(f14,'%d %d %f %f %f %f %f %f\n',[8,np-2]);
                    pipe_pair_nodes(pbpp+1:pbpp+np-2,1:2) = data(1:2,1:np-2)';
                    pipe_pair_values(pbpp+1:pbpp+np-2,1:6) = data(3:8,1:np-2)';
                end
                if(np > 1)
                    %Read line np
                    pos = ftell(f14);
                    data = fgets(f14);
                    if isempty(str2num(data))
                        fseek(f14,pos,'bof');
                        data = fscanf(f14,'%d %d %f %f %f %f %f %f\n',[8,1]);
                        fgets(f14); %Alpha numeric junk
                    else
                        %No alpha numeric junk present
                        data = str2num(data);
                    end
                    pipe_pair_nodes(pbpp+np-1,1:2) = data(1:2);
                    pipe_pair_values(pbpp+np-1,1:6) = data(3:8);
                end
                %Advance the counter
                pbpp = pbpp + np;
            else
                %Yet to be defined boundary types
            end
        end

        bndry_struct.name = [gridname,'_Boundary'];
        bndry_struct.nope = nelevbndy;
        bndry_struct.neta = nnodes_elev;
        bndry_struct.nvdllk = elev_seg;
        if (exist('elev_type','var')==1)
            bndry_struct.elev_type = elev_type;
        else
            bndry_struct.elev_type = [];
        end
        bndry_struct.elev_nodes = elev_nodes;
        bndry_struct.nbou = nflowbndy;
        bndry_struct.nvel = nnodes_flow;
        bndry_struct.nvellk = flow_seg;
        bndry_struct.flow_type = flow_type;
        bndry_struct.flow_nodes = flow_nodes;
        bndry_struct.weir_nodes = weir_nodes;
        bndry_struct.weir_values = weir_values;
        bndry_struct.weir_pair_nodes = weir_pair_nodes;
        bndry_struct.weir_pair_values = weir_pair_values;
        bndry_struct.pipe_pair_nodes = pipe_pair_nodes;
        bndry_struct.pipe_pair_values = pipe_pair_values;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp('No boundary information was found in the file.')
    bndry_struct = [];
end
%Close the file
fclose(f14);

return
