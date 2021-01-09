function obj = Make_f20_mycode(obj,filename,ts,te)

[edgesize,pertg,nvell_num] = Riverflux_distribution(obj);

data_riverflow = xlsread(filename);
river_volumeflow = data_riverflow(:,5:9);
time = data_riverflow(:,4);

total_node_num = sum(nvell_num);
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