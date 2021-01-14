function obj = Make_f20_volume_flow(obj,filename,ts,te,DT)
% obj = Make_f20_volume_flow(obj,filename,ts,te,DT)
% Input a msh class object to get the values of the nodal fluxes (q) (m^2/s) for 
% times between ts and te based on the input file. filename is the name of a column
% delimited file (csv) that stores the total volume flow (Q) (m^3/s) time series of  
% the cross-sections where the riverine boundaries in the mesh are located. 
%
% ts and te represent the start and end time of the total volume flow, respectively.
%
% The columns of the csv file should be organized as follows:  
% year,month,day,hour,minute,second,volume_flow_1,volume_flow_2... 
% The column order of the total volume flow (volume_flow_1,volume_flow_2...) 
% must be specified in the order in which the riverine boundaries appear in the 
% fort.14 file, or in which you make the riverine boundaries with the data cursor 
% method in msh.make_bc.
%
% A description of transforming the volume flux (m^3/s) to areal flux (m^2/s):
%
% The function of Riverflux_distribution (has been placed in the utilities folder)  
% will be called first to calculate the representative edge width and flux percentage  
% for each node on the riverine boundary. Nodal flux (m^3/s) can be calculated by
% multiplying the total volume flow (Q) (m^3/s) by the percentage of flow dedicated to 
% that boundary node. Then the nodal flux (q) (m^2/s) can be calculate by dividing 
% nodal flux by the nodal edge width.
%
% Note: the default total volume flow (Q) (m^3/s) is an hourly average time series, 
% thus the DT equals to 3600 (s). If your Q is an daily average time series, you should 
% use a value of 86400(24*3600) (s) for DT.
%
% Author:      Jiangchao Qiu                                
% Created:     January 7 2021                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[edgesize,pertg,nvell_num] = Riverflux_distribution(obj);

total_node_num = sum(nvell_num);
N = length(nvell_num);
 
data_riverflow = csvread(filename);

river_volumeflow = data_riverflow(:,7:7+N-1);

if DT==3600 % the total volume flow is an hourly average time series:
    time = data_riverflow(:,4);
elseif DT==86400 % the total volume flow is a daily average time series:
    time = data_riverflow(:,3);
else
    error('The total volume flow should be an hourly or daily average time series')
end

T = length(time);
F = zeros(total_node_num,T);

for i=1:T
    flow_this_time = river_volumeflow(i,:);
    flow_this_time1 = repmat(flow_this_time,[size(pertg,1),1]);
    flow_Q = (flow_this_time1.*pertg);
    flow_q = (flow_this_time1.*pertg)./edgesize;
    flow_q1 = flow_q(:);
    %for the condition that river boundaries have different number of nodes
    flow_q2= flow_q1(~isnan(flow_q1));  
    F(:,i) = flow_q2(:);  
end

%% Make into f20 struct
obj.f20.DataTitle = [datestr(ts) ' -> ' datestr(te)];
obj.f20.FTIMINC = DT;
obj.f20.NumOfNodes = size(F,1);
obj.f20.NumOfTimes = size(F,2);
obj.f20.Val = F(:);
%EOF
