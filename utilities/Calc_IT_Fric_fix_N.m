%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Title:       Calc_IT_Fric_fix_N                                 %
%  Description: The IT calculation process sometimes gets stuck    % 
%               in an infinite loop when using the original N data.%
%               This function is used to fill NaNs in the          %
%               originally generated N data in the vertical        %
%               direction by extrapolating all N values down the   %
%               column using the nearest non-NaN value above       %
%  Inputs:      Original Gridded_N_values.mat dataset              %            
%               (which is generated using Make_Gridded_N.m)        %     
%  Outputs:     An updated .mat file of profiles of N at each      %
%               desired contour level as specified at various      %
%               locations as the original Gridded_N_values.mat     %
%  Author:      Jiangchao Qiu                                      %
%  Created:     Feb 14 2023                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;
close all;  

%% load your original Gridded_N_values.mat datasey
load('H:\QJC\8-Github\OceanMesh2D_Jch\datasets\Gridded_N_values.mat')

% You can the choose to update the global or regional N values according 
% to your needs. My study area is in the Bay of Bangal, 
% the longitude is from 60 to 110,(so i is from 960 to 1160)
% the latitude is from 0 to 30,(so j is from 361 to 501）

% if you want to update N galobally, 
% i should be 1:1440 and j should be 1:720
for i = 960:1160
    for j = 361:501
        temp_vertical = N(i,j,:);
        temp_vertical_new = fillmissing(temp_vertical,'previous');
        N(i,j,:) = temp_vertical_new;
    end
end

o_name      = 'Gridded_N_values_new';

save([o_name '.mat'],'lon','lat','z','N');

