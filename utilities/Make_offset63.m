function obj = Make_offset63(obj,time_vector,offset_nodes,offset_values)
% obj = Make_offset63(obj,time_vector,offset_nodes,offset_values)
% 
% Inputs:   obj            - msh class 
%           time_vector    - datetime vector : size NTx1 or 1xNT
%           offset_nodes   - offset nodes    : size 1xNP
%           offset_values  - offset values   : size NTxNP
%
%  Author:      William Pringle                                 
%  Created:     Nov 22 2021                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(offset_nodes,1) ~= 1
   if size(offset_nodes,2) == 1 
       offset_nodes = offset_nodes';
   else
       error('offset_nodes size must be 1xNP')
   end
end
if size(offset_values,1) ~= length(time_vector)
   error('offset_values size of 1st dimension must be same as time_vector length')
end
if size(offset_values,2) ~= length(offset_nodes)
   error('offset_values size of 2nd dimension must be same as offset_nodes length')
end
        
% enter in msh.offset63 struct
obj.offset63.header = ['#dynamicwaterleveloffsets from ' ...
                datestr(time_vector(1)) ' to  ' datestr(time_vector(end))];
obj.offset63.num_times = length(time_vector);
obj.offset63.time_interval = round(...
                              seconds(time_vector(end) - time_vector(1))/ ...
                              (length(time_vector)-1)...
                             );
obj.offset63.default_val = 0;
obj.offset63.offset_nodes = offset_nodes;
obj.offset63.offset_values = offset_values;
%EOF
end


