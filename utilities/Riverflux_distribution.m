function [result_edgesize,result_pertg,nvell_num] = Riverflux_distribution(obj)     
% [result_edgesize,result_pertg，nvell_num] = Riverflux_distribution(obj) 
% 
% Input a msh class object to get the representative edge width and flux 
% percentage for each node on the riverine flow boundary.
% 
% This function is used to distribute a total flux of a cross-section where 
% the riverine boundary located to each node of the riverine boundary. 
%
% result_edgesize and result_pertgis are the result of the representative 
% edge width and flux percentage for each node, respectively. nvell_num 
% represents the number of nodes for each riverine boundary.
%
% Each column of result_edgesize and result_pertg represents the result of 
% each riverine boundary. 
%
% Nodal representative edge width equals to the sum of half the width of
% each of the two edges it is connected to.
%
% Nodeal flow percentage on the riverine boundary equals to the flow area
% of this node divided by the total flow area of the cross-section.
% 
% Nodal flow area is calculated via a trapezoidal rule using its representative
% edge width and the bathymetric depths of this node and its neighboring nodes.
%
% Note that the calculation of representative edge width and flow percentage 
% for the nodes at either end of the boundary is a little bit of different. 
% 
% User need to specify river flow boundaries (ibtype=22) before using it.
%
% Author:      Jiangchao Qiu                                
% Created:     January 7 2021                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bd_dat = obj.bd;
p_dat = obj.p;
b_dat = obj.b;
river_num = find(bd_dat.ibtype == 22);
N = length(river_num);% total number of river boundarys
nvell_num = bd_dat.nvell(river_num);% number of nodes for each river boundary

% Consider the riverine boundaries may have different number of nodes, the 
% number of rows for result_edgesize and result_pertg is set to max(nvell_num)
result_edgesize = NaN(max(nvell_num),N); % the final result of edgesize
result_pertg = NaN(max(nvell_num),N);% the final result of percentage

if isempty(river_num)
    error('No riverine boundary information to distribute total flow')
end

for i=1:N
    node_num = bd_dat.nbvv(:,river_num(i));
    node_num(node_num==0)=[];
    bathy = b_dat(node_num);
    location = p_dat(node_num,:); 
    J = length(node_num);% total number of nodes for the current river boundary
    flowarea = zeros(J,1);
    edgesize = zeros(J,1);
    %% calculate the projected distance for each edge on the boundary
    proj = 'Mercator';
    m_proj(proj,'lon',[ min(obj.p(:,1))-0.25 max(obj.p(:,1))+0.25 ],...
                'lat',[ min(obj.p(:,2))-0.25 max(obj.p(:,2))+0.25])
    proj_location = zeros(J,2);
    for j=1:J
        [proj_location(j,1),proj_location(j,2)] = m_ll2xy(location(j,1),location(j,2));
    end
    node_distance = 1000*m_xydist(proj_location(:,1),proj_location(:,2));     
    %% calculate the respresentive width for each node
    % note: if the current boundary has J nodes, there are J-1 edges	
    for j=1:J
		if j==1      % the first edge (the number is J) 
            edgesize(j) = 0.5*node_distance(j);
        elseif j==J  % the last edge (the number is J-1)
            edgesize(j) = 0.5*node_distance(j-1);
        else         % the other edges 
            edgesize(j) = 0.5*node_distance(j-1)+0.5*node_distance(j);
        end
        result_edgesize(1:J,i) = edgesize;
    end
    %% calculate flow area for each node on the boundary
    for j=1:J
        if j==1      % the first node
            local_z1 = (bathy(j)+bathy(j+1))*0.5;
            flowarea(j) = 0.5*(bathy(j)+local_z1)*(0.5*node_distance(j));
        elseif j==J  % the last node 
            local_z1 = (bathy(j-1)+bathy(j))*0.5;
            flowarea(j) = 0.5*(bathy(j)+local_z1)*(0.5*node_distance(j-1));
        else         % the other nodes
            local_z1 = (bathy(j-1)+bathy(j))*0.5;
            local_z2 = (bathy(j)+bathy(j+1))*0.5;
            flowarea(j) = 0.5*(bathy(j)+local_z1)*(0.5*node_distance(j-1))+...
                          0.5*(bathy(j)+local_z2)*(0.5*node_distance(j));
        end
        percent = flowarea/sum(flowarea);  % the flux distribution percentage for each node
        result_pertg(1:J,i) =  percent;
    end
end


