function obj = Make_f20_volumeflow(obj,filename,ts,te)
% obj = Make_f20_volumeflow(obj,filename,ts,te)
% Input a msh class object to get the values of the normal fluxes for times 
% between ts and te based on the input files. filename is the name of a csv 
% file that stores the volume flow (m3/s) data time series. the csv file should  
% be ut into the datasets folder. The first six column record the Year, Month, 
% Day, Hour, Minute, Second of the volume flow. The Station sequencing is assumed 
% to match what is specified in the riverine boundary setting. 
%
%  Author:      Jiangchao Qiu                                
%  Created:     January 7 2021                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[edgesize,pertg,nvell_num] = Riverflux_distribution(obj);

total_node_num = sum(nvell_num);
N = length(nvell_num);
 
data_riverflow = csvread(filename);

river_volumeflow = data_riverflow(:,7:7+N-1);
time = data_riverflow(:,4);

T = length(time);
DT = 3600; %[s] 
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