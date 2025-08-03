%% SIngle Subject COBE HCP + NKI use old dictionaries
clear
working_dir = 'H:\Single_scan_cobe_results2';
file_path_nki = fullfile(working_dir,'NKI_matfiles'); %NKI average timeseries path

addpath(fullfile(working_dir,'codes_functions','Dynamic_fingerprints-main')) % load amico dfc idiff codes
addpath(fullfile(working_dir,'codes_functions','functions'))
addpath(fullfile(working_dir,'codes_functions','functions','COBE'))
%% Load the matfiles

net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks

matfile_path = fullfile(working_dir,'matfiles'); 
atlas_name = 's400_hcp';
para_check_path = fullfile(working_dir,'para_check');
results_path = fullfile(working_dir,'from_hcp_dict_nki',atlas_name);
mkdir(results_path)

win_hcp = 400;
win_nki = 288;

st_hcp = 20;
st_nki = 7;

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

hcp_perm_matfiles_path = ['E:\dy_identification_cobe\new_data_results\perm_check\' atlas_name];
for i_net = 9:-1:1
    load(fullfile(hcp_perm_matfiles_path,[net_names{i_net},'_max_acco1_' atlas_name ,'_test5x2.mat']),'cobe_dict')


    n_nki = size(data_nki_1,3);
    nki_sub = randperm(n_nki);

    %% Get yeo network nodes
    % 1 --> node belongs to Visual network
    % 2 --> node belongs to Somatomotor network
    % 3 --> node belongs to Dorsal Attention network
    % 4 --> node belongs to Ventral Attention network
    % 5 --> node belongs to Limbic network
    % 6 --> node belongs to Fronto-Parietal network
    % 7 --> node belongs to Default Mode network
    % 0 --> node does not belongs to any of the 7 yeo network

    atlas_path = fullfile(working_dir,'atlases');
    atlas_brain_path = fullfile(atlas_path,[atlas_name,'.nii']);

    atlas_yeo_path = fullfile(atlas_path,'yeo_7_net.nii');
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


    disp([num2str(i_net) 'network'])
    if i_net == 9
        clear net_id
    else
        net_id = nets(i_net);
    end

    if exist('net_id','var')                                                    % if net_id exists
        net = idx == net_id;                                                    % Get the index of the ROIs corresponding to the Yeo network
    else
        net = 1:size(data_nki_1,2);                                                 % else use all the ROIs
    end

    configs.Nparc = nnz(net);                                                   % Number of brain regions
    configs.mask_ut = triu(true(configs.Nparc,configs.Nparc),1);                % Upper triangular mask

    run1_fs = filter(b_nki,a_nki,data_nki_1(:,net,nki_sub));                                   % filter the run1_lr data
    vec_nki_1 = functional_connectivity(run1_fs,0,win_nki,st_nki);               % DFC computation for run 1

    run2_fs = filter(b_nki,a_nki,data_nki_2(:,net,nki_sub));                                   % filter the run1_lr data
    vec_nki_2 = functional_connectivity(run2_fs,0,win_nki,st_nki);               % DFC computation for run 1



    vec_test_nki_1 = vec_nki_1;
    vec_test_nki_2 = vec_nki_2;

    configs.k_frames = 1:size(vec_test_nki_1,2);                          % number of dfc frames
    for ii = 1:5
        for cv = 1:2

        %test
        Sub_spec_test_1_nki = cobe_transform(vec_test_nki_1,cobe_dict(:,:,cv,ii));                                    % the dictionary extracted from the training subjects.
        Sub_spec_test_2_nki = cobe_transform(vec_test_nki_2,cobe_dict(:,:,cv,ii));
        s = size(vec_test_nki_1,2);                                           % number of frames in dFC
        p = size(vec_test_nki_1,3);
        % Original dI_diff
        configs.k_frames = 1:size(vec_test_nki_1,2);                          % number of dfc frames
        [dID_test_nki,k_Frames_run1_test_nki,k_Frames_run2_test_nki,dIdent_test_nki] = f_compute_Idiff_frames(permute(vec_test_nki_1,[3,1,2]),permute(vec_test_nki_2,[3,1,2]),configs);
        [~,max_dI_idx_test_nki] = max(dIdent_test_nki);

        [Idiff_max_test_perm_nki(cv,ii),Idiff_max_test_idx_perm_nki] = max(dID_test_nki.Idiff);

        % Subject-specific dI_diff
        [dID_cobe_test_nki,k_Frames_run1_cobe_test_nki,k_Frames_run2_cobe_test_nki,dIdent_cobe_test_nki] = f_compute_Idiff_frames(permute(Sub_spec_test_1_nki,[3,1,2]),permute(Sub_spec_test_2_nki,[3,1,2]),configs);
        [~,max_dI_idx_cobe_test_nki] = max(dIdent_cobe_test_nki);

        [Idiff_max_test_cobe_perm_nki(cv,ii),Idiff_max_test_cobe_idx_perm_nki] = max(dID_cobe_test_nki.Idiff);
        %     Classification rate all
        [acco_test_perm_nki(cv,ii),acco_or_test_perm_nki(cv,ii)] = classification_rate_dfc(max_dI_idx_cobe_test_nki',max_dI_idx_test_nki',size(vec_test_nki_1,3),size(vec_test_nki_1,2));


        % figure; fc_imshow(veccorr(cobe_dict(:,1)),sort_idx_val,net_names)

        disp('#################################################################################################')
       
        disp('---------------------------------------------------------------------------')
        disp(['Original dIdiff (nki): ' num2str(Idiff_max_test_perm_nki(cv,ii))])
        disp(['COBE dIdiff (nki): ' num2str(Idiff_max_test_cobe_perm_nki(cv,ii))])

     
        disp('-------------------------------------------------------------------------------------------------')
        disp(['Original acco (nki): ' num2str(acco_or_test_perm_nki(cv,ii))])
        disp(['COBE acco (nki): ' num2str(acco_test_perm_nki(cv,ii))])

        disp('#################################################################################################')

        save(fullfile(results_path,[char(net_names(i_net)),'_nki_' atlas_name '_test_from_hcp_perm_dict_5x2.mat']),"Idiff_max_test_perm_nki",...  _cobe_perm_acco_idiff_weight_',num2str(w_idiff*100),'_max2.mat
            "Idiff_max_test_cobe_perm_nki",'acco_or_test_perm_nki','acco_test_perm_nki',"cobe_dict",'win_nki','st_nki','-v7.3')
        end
    end
end

function [acco,acco_or] = classification_rate_dfc(idf_idx,idf_idx_orr,p,s)

%% Input
% sub_spec --> subject specific FC (train) --> n x p (n = FC upper triangular vector, p = number of training scans)
% sub_spec_test --> subject specific FC (test) --> n x p (p = number of testing scans )
% maph_n --> normalized original FC matrix used to train the dictionary --> n x p1
% maph1_n --> normalized original FC matrix used for testing --> n x p2
% s = sessions per subjects
gt = repmat(1:p,s,1);   % Ground truth
gt = reshape(gt,[],1);   % Ground truth

out_or = reshape(gt(idf_idx_orr),s,[]);
out_mode_or = mode(out_or);

out_ = reshape(gt(idf_idx),s,[]);
out_mode_ = mode(out_);
gt = 1:p;
confu_ = confusionmat(gt,out_mode_);
confu_orr = confusionmat(gt,out_mode_or);

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

