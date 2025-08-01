%% Get Average Parameters value
% This script was written by Pratik Jain if you have any questions email
% jainpratik412@gmail.com
% 
% This Script uses the results of parameter check from previous script
% parameter_search.m and averages the optimal parameters across the blocks
clear

working_dir = 'H:\Single_scan_cobe_results2';


%% Load the matfiles
matfile_path = 'E:\dy_identification_cobe\matfiles';                         % location of the matfiles containg the region timeseries extracted from atlas

results_path = fullfile(working_dir,'para_check');
atlas_name = 's400_hcp';

Win_size = [50,100,200,400,800];
Com_comp = 1:2:20;
Stride = [20,50,100];
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'};



for i_net = 1:9
    acco_avg = zeros(length(Win_size),length(Stride),length(Com_comp));
    acco_avg_or = zeros(length(Win_size),length(Stride),length(Com_comp));
    for i_block = 1:3
        
        load(fullfile(results_path,[atlas_name,'_para_test',num2str(i_block),'_',char(net_names{i_net}),'.mat']))
        
        acco_avg = acco_avg + acco;

        acco_avg_or = acco_avg_or + acco_or;

    end
   
    acco_avg = acco_avg./3;    
    acco_avg_or = acco_avg_or./3; 
    save(fullfile(results_path,[atlas_name,'_para_test_avg',char(net_names{i_net}),'.mat']),'acco_avg_or','acco_avg');
end
