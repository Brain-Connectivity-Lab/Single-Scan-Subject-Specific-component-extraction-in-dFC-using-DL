%% SIngle Subject COBE HCP + NKI
clear
% rng(5,"twister");
%% Load the matfiles
addpath('E:\dy_identification_cobe\Dynamic_fingerprints-main') % load amico dfc idiff codes
addpath('D:\codes\Individual_differences_code_HBM_paper\Codes\functions')
addpath('C:\Users\krish\Box\Single Subject COBE\paper_figure_codes_final\functions\COBE')
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks

matfile_path = 'E:\dy_identification_cobe\matfiles';%'E:\hcp_100_unrelated_matfiles\removed_nan_subject_from_all';                         % location of the matfiles containg the region timeseries extracted from atlas
atlas_name = 's400_hcp';
para_check_path = 'E:\dy_identification_cobe\new_data_results\para_check';
results_path = ['C:\Users\krish\Box\Single Subject COBE\new_data_results\nki_hcp\' atlas_name '\seperate'];
mkdir(results_path)
file_path_nki = 'D:\DATA\Mat_files\NKI\';

n_train_hcp = 41;
n_train_nki = 41;

% Win_size = [50,100,200,400,800];
% Com_comp = 1:2:20;
% Stride = [20,50,100];
win_hcp = 400;
win_nki = 206;

st_hcp = 20;
st_nki = 10;

Component_remv = 5;
%% HCP block 1
% load(fullfile('E:\dy_identification_cobe\parameter_check\',"s400_hcp_new_my_cobec.mat"))                            % parameter check file
data_hcp_1 = [];
data_hcp_2 = [];
subname_1 = cell(1);
subname_2 = cell(1);

for i_block = 1
    %% Load Rest data

    data_1 = load(fullfile(matfile_path,['HCP_' atlas_name '_run_1_encode_LR_sub_block' num2str(i_block) '_reg_ts.mat']),'reg_ts','p_sub_name');  % load the run1 LR subjects
    sub_namerun1 = data_1.p_sub_name;                                                                                  % run1 LR subject names
    data_1 = data_1.reg_ts;                                                                                            % run1 LR average region time series


    data_2 = load(fullfile(matfile_path,['HCP_' atlas_name '_run_1_encode_RL_sub_block' num2str(i_block) '_reg_ts.mat']),'reg_ts','p_sub_name');  % load the run1 LR subjects
    sub_namerun2 = data_2.p_sub_name;                                                                                  % run1 LR subject names
    data_2 = data_2.reg_ts;                                                                                               % run2 LR average region time series


    data_hcp_1 = cat(3,data_hcp_1,data_1);
    data_hcp_2 = cat(3,data_hcp_2,data_2);

    subname_1 = cat(1,subname_1,sub_namerun1);
    subname_2 = cat(1,subname_2,sub_namerun2);

    clear data_1 data_2

end
subname_1(1) = [];
subname_2(1) = [];



configs.numSubj = 41;%(2*n_train_hcp - n_train_hcp) + (n_nki - n_train_nki);

%% NKI 94 subjects

data_1 = load([file_path_nki,'session_1_reg_ts_' atlas_name '_rest'],'reg_ts','p_sub_name'); % load the FLU1 subjects
sub_namedata1 = data_1.p_sub_name;                                              % FLU1 subject names
data_nki_1 = data_1.reg_ts;                                                    % FLU1 average region time series

data_2 = load([file_path_nki,'session_2_reg_ts_' atlas_name '_rest'],'reg_ts','p_sub_name');             % load the BAS1 subjects
sub_namedata2 = data_2.p_sub_name;                                              % BAS1 subject names
data_nki_2 = data_2.reg_ts;                                                    % BAS1 average region time series

sub_info_rest = cell(10,1);                                                      % will contain selected subject name and its age

[nki_sub_name,ia,ib] = intersect(string(sub_namedata1(:,1)),string(sub_namedata2(:,1)));

data_nki_1 = data_nki_1(:,:,ia);
data_nki_2 = data_nki_2(:,:,ib);


for i_net = [1,2,7,3,4,5,6,8,9]
    clear cobe_dict
    for ii = 1:20
        n_hcp = size(data_hcp_1,3);
        n_nki = size(data_nki_1,3);
        hcp_sub = randperm(n_hcp,2*n_train_hcp);
        nki_sub = randperm(n_nki,2*n_train_nki);

        %% Get yeo network nodes
        % 1 --> node belongs to Visual network
        % 2 --> node belongs to Somatomotor network
        % 3 --> node belongs to Dorsal Attention network
        % 4 --> node belongs to Ventral Attention network
        % 5 --> node belongs to Limbic network
        % 6 --> node belongs to Fronto-Parietal network
        % 7 --> node belongs to Default Mode network
        % 0 --> node does not belongs to any of the 7 yeo network

        atlas_brain_path = ['D:\Atlases\2mm_new\' atlas_name '.nii'];                      % Path of brain atlas used to get the Averaged time series
        atlas_yeo_path = 'D:\Atlases\2mm_new\yeo_7_net.nii';
        idx = yeo_networks(atlas_brain_path,atlas_yeo_path);                       % Identify the networks indexes for each ROI
        [sort_idx_val,sort_idx] = sort(idx);
        % net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks
        nets = 0:7;       % Choose yeo network, if net_id does not exist code will consider all the nodes (whole brain)

        col = [0 1 1
            0.4940 0.1840 0.5560;
            0 0.4470 0.4710;
            0 1 0;
            0 0 0;
            0.9290 0.694 0.125;
            1 0 1;
            1 0 0];
        %% Filter
        fc1 = 0.01;                                                                 % Low pass cutoff
        fc2 = 0.1;                                                                  % High pass cutoff
        TR_nki = 1.4;                                                                   % TR
        fs_nki = 1/TR_nki;                                                                  % fs
        TR_hcp = 0.72;                                                                   % TR
        fs_hcp = 1/TR_hcp;                                                                  % fs

        [b_nki,a_nki] = butter(5,[fc1,fc2]/(fs_nki/2),'bandpass');                              % butterworth filter of the 5th order
        [b_hcp,a_hcp] = butter(5,[fc1,fc2]/(fs_hcp/2),'bandpass');                              % butterworth filter of the 5th order


        disp([num2str(i_net) 'network'])
        if i_net == 9
            clear net_id
        else
            net_id = nets(i_net);
        end

        if exist('net_id','var')                                                    % if net_id exists
            net = idx == net_id;                                                    % Get the index of the ROIs corresponding to the Yeo network
        else
            net = 1:size(data_hcp_1,2);                                                 % else use all the ROIs
        end

        configs.Nparc = nnz(net);                                                   % Number of brain regions
        configs.mask_ut = triu(true(configs.Nparc,configs.Nparc),1);                % Upper triangular mask

        run1_fs = filter(b_hcp,a_hcp,data_hcp_1(:,net,hcp_sub));                                   % filter the run1_lr data
        vec_hcp_1 = functional_connectivity(run1_fs,0,win_hcp,st_hcp);               % DFC computation for run 1

        run2_fs = filter(b_hcp,a_hcp,data_hcp_2(:,net,hcp_sub));                                   % filter the run1_lr data
        vec_hcp_2 = functional_connectivity(run2_fs,0,win_hcp,st_hcp);               % DFC computation for run 1

        run1_fs = filter(b_nki,a_nki,data_nki_1(:,net,nki_sub));                                   % filter the run1_lr data
        vec_nki_1 = functional_connectivity(run1_fs,0,win_nki,st_nki);               % DFC computation for run 1

        run2_fs = filter(b_nki,a_nki,data_nki_2(:,net,nki_sub));                                   % filter the run1_lr data
        vec_nki_2 = functional_connectivity(run2_fs,0,win_nki,st_nki);               % DFC computation for run 1


%         vec_train = cat(3,vec_hcp_1(:,:,1:n_train_hcp),vec_nki_1(:,:,1:n_train_nki));
        vec_train_cell = cell(1);

        for i_n = 1:size(vec_hcp_1,3)
            vec_train_cell{i_n} = vec_hcp_1(:,:,i_n);
        end
%         for i_n = n_train_hcp+1:n_train_hcp + n_train_nki
%             vec_train_cell{i_n} = vec_nki_1(:,:,i_n-n_train_hcp);
%         end

        vec_test_nki_1 = vec_nki_1;
        vec_test_nki_2 = vec_nki_2;


        %% COBE

        COBE_CM_Component = my_cobec_cell( vec_train_cell,Component_remv ); % COBE algorithm



        Sub_spec_test_1_nki = cobe_transform(vec_test_nki_1,COBE_CM_Component);                                    % the dictionary extracted from the training subjects.
        Sub_spec_test_2_nki = cobe_transform(vec_test_nki_2,COBE_CM_Component);

       
        % Original dI_diff
        configs.k_frames = 1:size(vec_test_nki_1,2);                          % number of dfc frames
        [dID_test_nki,k_Frames_run1_test_nki,k_Frames_run2_test_nki,dIdent_test_nki] = f_compute_Idiff_frames(permute(vec_test_nki_1,[3,1,2]),permute(vec_test_nki_2,[3,1,2]),configs);
        [~,max_dI_idx_test_nki] = max(dIdent_test_nki);

        [Idiff_max_test_perm_nki(ii),Idiff_max_test_idx_perm_nki] = max(dID_test_nki.Idiff);

        % Subject-specific dI_diff
               
        configs.k_frames = 1:size(vec_test_nki_1,2);                          % number of dfc frames
        [dID_cobe_test_nki,k_Frames_run1_cobe_test_nki,k_Frames_run2_cobe_test_nki,dIdent_cobe_test_nki] = f_compute_Idiff_frames(permute(Sub_spec_test_1_nki,[3,1,2]),permute(Sub_spec_test_2_nki,[3,1,2]),configs);
        [~,max_dI_idx_cobe_test_nki] = max(dIdent_cobe_test_nki);

        [Idiff_max_test_cobe_perm_nki(ii),Idiff_max_test_cobe_idx_perm_nki] = max(dID_cobe_test_nki.Idiff);
        %     Classification rate all
        [acco_test_perm_nki(ii),acco_or_test_perm_nki(ii)] = classification_rate_dfc(max_dI_idx_cobe_test_nki',max_dI_idx_test_nki',size(vec_test_nki_1,3),size(vec_test_nki_1,2));


        cobe_dict(:,:,ii) = COBE_CM_Component;


        % figure; fc_imshow(veccorr(cobe_dict(:,1)),sort_idx_val,net_names)


        disp('---------------------------------------------------------------------------')
        disp(['Original dIdiff (nki): ' num2str(Idiff_max_test_perm_nki(ii))])
        disp(['COBE dIdiff (nki): ' num2str(Idiff_max_test_cobe_perm_nki(ii))])

        disp('-------------------------------------------------------------------------------------------------')
        disp(['Original acco (nki): ' num2str(acco_or_test_perm_nki(ii))])
        disp(['COBE acco (nki): ' num2str(acco_test_perm_nki(ii))])

        disp('#################################################################################################')

        save(fullfile(results_path,[char(net_names(i_net)),'_nki_hcp_seperate' atlas_name '_test.mat']),"Idiff_max_test_perm_nki",...  _cobe_perm_acco_idiff_weight_',num2str(w_idiff*100),'_max2.mat
            "Idiff_max_test_cobe_perm_nki",'acco_or_test_perm_nki','acco_test_perm_nki',"cobe_dict",'Component_remv','win_hcp','st_hcp','win_hcp','st_hcp','-v7.3')
    end
end
function [acco,acco_or] = classification_rate_dfc(idf_idx,idf_idx_orr,p,s)

% Input
% % sub_spec --> subject specific FC (train) --> n x p (n = FC upper triangular vector, p = number of training scans)
% % sub_spec_test --> subject specific FC (test) --> n x p (p = number of testing scans )
% % maph_n --> normalized original FC matrix used to train the dictionary --> n x p1
% % maph1_n --> normalized original FC matrix used for testing --> n x p2
% s = sessions per subjects
gt = repmat(1:p,s,1);   % Ground truth
gt = reshape(gt,[],1);   % Ground truth

confu_ = confusionmat(gt,gt(idf_idx));
confu_orr = confusionmat(gt,gt(idf_idx_orr));

acco = sum(diag(confu_))/sum(confu_(:));
acco_or = sum(diag(confu_orr))/sum(confu_orr(:));
end

function [Sub_spec_test] = cobe_transform(vec,COBE_CM_Component)

flat2d_run = reshape(vec,size(vec,1),size(vec,2)*size(vec,3));
X_run_test = COBE_CM_Component'*flat2d_run;
% flat2d_Com_run_test = COBE_CM_Component*X_run_test;
flat2d_Sub_spec_run_test = flat2d_run - COBE_CM_Component*X_run_test;

Sub_spec_test = reshape(flat2d_Sub_spec_run_test,size(vec));
% Comm_test = reshape(flat2d_Com_run_test,size(vec));
end

