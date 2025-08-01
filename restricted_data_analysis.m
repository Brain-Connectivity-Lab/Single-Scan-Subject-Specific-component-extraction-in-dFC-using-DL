%% Twin and family data
clear
working_dir = 'H:\Single_scan_cobe_results2';

addpath(fullfile(working_dir,'codes_functions','Dynamic_fingerprints-main')) % load amico dfc idiff codes
addpath(fullfile(working_dir,'codes_functions','functions'))
addpath(fullfile(working_dir,'codes_functions','functions','COBE'))

matfile_path = fullfile(working_dir,'matfiles'); 
atlas_name = 's400_hcp';
%% Atlas
res_type = 'MZ';
para_check_path = fullfile(working_dir,'para_check');
results_path = fullfile(working_dir,'5x2cv_twins',atlas_name,res_type);
mkdir(results_path)
Win_size = [50,100,200,400,800];
Com_comp = 1:2:20;
Stride = [20,50,100];
at_names = {'a116','d164','p264','w300','s400_hcp'}; %AAL, Dosenbach, Power, Seitzman, Shaefer_HCP 
run = '1';
encode= 'RL';
% Load all data
reg_ts1 = [];
p_sub_name1 = [];
for i = 1:3
    data = load(fullfile(matfile_path,['HCP_',atlas_name,'_run_',run,'_encode_',encode,'_sub_block',num2str(i),'_reg_ts.mat']));
    p_sub_name1 = [p_sub_name1;data.p_sub_name];
    reg_ts1 = cat(3,reg_ts1,data.reg_ts);
end
clear data

encode = 'LR';
reg_ts2 = [];
p_sub_name2 = [];
for i = 1:3
    data = load(fullfile(matfile_path,['HCP_',atlas_name,'_run_',run,'_encode_',encode,'_sub_block',num2str(i),'_reg_ts.mat']));
    p_sub_name2 = [p_sub_name2;data.p_sub_name];
    reg_ts2 = cat(3,reg_ts2,data.reg_ts);
end
clear data

%% Restricted csv

csv = readtable('C:\Users\krish\Box\Single Subject COBE\RESTRICTED_pratikpj44_3_21_2025_15_43_25_new.csv');

[~, ia] = unique(csv.Family_ID);
% 
unrelated_fam = csv(ia,:);
[~, ia] = unique(unrelated_fam.Mother_ID);
unrelated_fam_m = unrelated_fam(ia,:);
[unrelated_, ia] = unique(unrelated_fam_m.Father_ID);
unrelated_fam_m_f = unrelated_fam_m(ia,:);

twin_mz = csv(string(csv.Zygosity_preffered) == 'MZ',:); %Zygosity_preffered
twin_dz = csv(string(csv.Zygosity_preffered) == 'DZ',:); %Zygosity_preffered
switch res_type
    case 'MZ'
        final_csv = twin_mz;
        
    case 'un'
        final_csv = unrelated_fam_m_f;

    case 'nt_un'
        final_csv = unrelated_fam_m_f;
        final_csv = final_csv(string(table2cell(final_csv(:,4)))=='NotTwin',:);
    case 'DZ'
        final_csv = twin_dz;

end
%% Intersect CSV and reg_ts

[~,ia,ib] = intersect(string(p_sub_name1),string(table2cell(final_csv(:,1))));

reg_ts1 = reg_ts1(:,:,ia);
reg_ts2 = reg_ts2(:,:,ia);
p_sub_name = p_sub_name1(ia);
final_csv = final_csv(ib,:);

switch res_type
    case 'MZ'
        TBL = tabulate(final_csv.Family_ID);
        TBL(cell2mat(TBL(:,2))~= 1,:) = [];
        [temp,ia2,ib2] = intersect(string(TBL(:,1)),string(table2cell(final_csv(:,7))));
        final_csv(ib2,:) = [];
        reg_ts1(:,:,ib2) = [];
        reg_ts2(:,:,ib2) = [];
        p_sub_name(ib2) = [];

        [~,idf] = sort(final_csv.Family_ID);
        final_csv = final_csv(idf,:);
        reg_ts1 = reg_ts1(:,:,idf);
        reg_ts2 = reg_ts2(:,:,idf);
        p_sub_name = p_sub_name(idf);
    case 'DZ'
        TBL = tabulate(final_csv.Family_ID);
        TBL(cell2mat(TBL(:,2))~= 1,:) = [];
        [temp,ia2,ib2] = intersect(string(TBL(:,1)),string(table2cell(final_csv(:,7))));
        final_csv(ib2,:) = [];
        reg_ts1(:,:,ib2) = [];
        reg_ts2(:,:,ib2) = [];
        p_sub_name(ib2) = [];

        [~,idf] = sort(final_csv.Family_ID);
        final_csv = final_csv(idf,:);
        reg_ts1 = reg_ts1(:,:,idf);
        reg_ts2 = reg_ts2(:,:,idf);
        p_sub_name = p_sub_name(idf);

end

confu_orr = [];
confu_orr2 = [];
confu_  = [];
confu_2 = [];
% will contain selected subject name
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

net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks
nets = 0:7;       % Choose yeo network, if net_id does not exist code will consider all the nodes (whole brain)

col = [0 1 1
    0.4940 0.1840 0.5560;
    0 0.4470 0.4710;
    0 1 0;
    0 0 0;
    0.9290 0.694 0.125;
    1 0 1;
    1 0 0];


for i_net =7%9:-1:1
    % during train and test for comparision
    clear cobe_dict
    for iter = 1:5  % 5 iterations
        rng(iter,"twister");
        rand_sess = randperm(2,2);
        if strcmp(res_type,'MZ') || strcmp(res_type,'DZ')
            k_fold_cv = cvpartition(length(p_sub_name)/2,"KFold",2);
        else
            k_fold_cv = cvpartition(length(p_sub_name),"KFold",2);
        end
        for i = 1:2
            if strcmp(res_type,'MZ') || strcmp(res_type,'DZ')
                tt1 = test(k_fold_cv,i);
                tt = [tt1';tt1'];
                idxcv(:,i) = tt(:);
                rand_sess_train(i) = rand_sess(1);
            else
                idxcv(:,i) = test(k_fold_cv,i);
                rand_sess_train(i) = rand_sess(1);
            end
        end

        %     close all
        disp([num2str(i_net) 'network'])

        para_check = load(fullfile(para_check_path,[atlas_name '_para_test_avg' char(net_names(i_net)) '.mat']));

        [val,idx1] = max(squeeze(para_check.acco_avg),[],'all');                                                               % Get the idx of the max training idiff value
        [win_i,S_i,C_i] = ind2sub(size(squeeze(para_check.acco_avg)),idx1);
        %     win_i = 3;
        Component_remv  = Com_comp(C_i);
        win = Win_size(win_i);%400;                                                                                           % window
        st = Stride(S_i);%50;

        if i_net == 9
            clear net_id
        else
            net_id = nets(i_net);
        end

        if exist('net_id','var')                                                    % if net_id exists
            net = idx == net_id;                                                    % Get the index of the ROIs corresponding to the Yeo network
        else
            net = 1:size(reg_ts1,2);                                                 % else use all the ROIs
        end

        configs.Nparc = nnz(net);                                                   % Number of brain regions
        configs.mask_ut = triu(true(configs.Nparc,configs.Nparc),1);                % Upper triangular mask


        %% Filter

        fc1 = 0.01;                                                                 % Low pass cutoff
        fc2 = 0.1;                                                                  % High pass cutoff
        TR = 0.72;                                                                   % TR
        fs = 1/TR;                                                                  % fs

        [b,a] = butter(5,[fc1,fc2]/(fs/2),'bandpass');                              % butterworth filter of the 5th order

%         cobe_dict = zeros(nnz(configs.mask_ut),Component_remv,2,5);


        for p_i = 1:2
            clear vec_rest
            
            if rand_sess_train(p_i) == 1
                run1_fs = filter(b,a,reg_ts1(:,net,idxcv(:,p_i)));                                   % filter the run1_lr data
                vec_rest(:,:,:) = functional_connectivity(run1_fs,0,win,st);               % DFC computation for run 1
            else
                run2_fs = filter(b,a,reg_ts2(:,net,idxcv(:,p_i)));                                   % filter the run2_lr data
                vec_rest(:,:,:) = functional_connectivity(run2_fs,0,win,st);               % DFC computation for run 2
            end

            if p_i == 1
                f_i = 2;
            elseif p_i == 2
                f_i = 1;
            end
            configs.numSubj = nnz(idxcv(:,f_i));
            COBE_CM_Component = my_cobec( vec_rest,Component_remv ); % COBE algorithm

            run1_fs = filter(b,a,reg_ts1(:,net,idxcv(:,f_i)));                                   % filter the run1_lr data
            clear vec_rest_test_1 vec_rest_test_2
            vec_rest_test_1(:,:,:) = functional_connectivity(run1_fs,0,win,st);               % DFC computation for run 1

            run2_fs = filter(b,a,reg_ts2(:,net,idxcv(:,f_i)));                                   % filter the run1_rl data
            vec_rest_test_2(:,:,:) = functional_connectivity(run2_fs,0,win,st);               % DFC computation for run 1
            configs.k_frames = 1:size(vec_rest_test_2,2);                          % number of dfc frames

            %% CobeIndiv_cell_run

            [Comm_test_1,Sub_spec_test_1] = cobe_transform(vec_rest_test_1,COBE_CM_Component);                                    % the dictionary extracted from the training subjects.
            [Comm_test_2,Sub_spec_test_2] = cobe_transform(vec_rest_test_2,COBE_CM_Component);

            % I_diff value after individual differentiability is extracted using COBE

            s = size(vec_rest_test_1,2);                                           % number of frames in dFC
            p = size(vec_rest_test_1,3);                                           % number of test subjects
            %% Original dI_diff
            [dID_test,k_Frames_run1_test,k_Frames_run2_test,dIdent_test] = f_compute_Idiff_frames(permute(vec_rest_test_1,[3,1,2]),permute(vec_rest_test_2,[3,1,2]),configs);
            [~,max_dI_idx_test] = max(dIdent_test);

            [Idiff_max_test_perm(p_i,iter),Idiff_max_test_idx_perm(p_i,iter)] = max(dID_test.Idiff);

            %% Subject-specific dI_diff
            [dID_cobe_test,k_Frames_run1_cobe_test,k_Frames_run2_cobe_test,dIdent_cobe_test] = f_compute_Idiff_frames(permute(Sub_spec_test_1,[3,1,2]),permute(Sub_spec_test_2,[3,1,2]),configs);
            [~,max_dI_idx_cobe_test] = max(dIdent_cobe_test);

            [Idiff_max_test_cobe_perm(p_i,iter),Idiff_max_test_cobe_idx_perm(p_i,iter)] = max(dID_cobe_test.Idiff);

            %% common

            [dID_cobe_test_com,k_Frames_run1_cobe_test_com,k_Frames_run2_cobe_test_com,dIdent_cobe_test_com] = f_compute_Idiff_frames(permute(Comm_test_1,[3,1,2]),permute(Comm_test_2,[3,1,2]),configs);
            [~,max_dI_idx_cobe_test_com] = max(dIdent_cobe_test_com);

            [Idiff_max_test_cobe_perm_com(p_i,iter),Idiff_max_test_cobe_idx_perm_com(p_i,iter)] = max(dID_cobe_test_com.Idiff);

            %%     Classification rate all
            if p_i == 1
                [acco_test_perm(p_i,iter),acco_or_test_perm(p_i,iter),confu_(:,:,p_i,iter),confu_orr(:,:,p_i,iter)] = classification_rate_dfc(max_dI_idx_cobe_test',max_dI_idx_test',p,s);
            elseif p_i == 2
                [acco_test_perm(p_i,iter),acco_or_test_perm(p_i,iter),confu_2(:,:,p_i,iter),confu_orr2(:,:,p_i,iter)] = classification_rate_dfc(max_dI_idx_cobe_test',max_dI_idx_test',p,s);
            end
            cobe_dict(:,:,p_i,iter) = COBE_CM_Component;

            disp('#################################################################################################')
            disp(i_net)
            disp(p_i)
            disp(iter)
            disp('---------------------------------------------------------------------------')
            disp(['Original dIdiff: ' num2str(Idiff_max_test_perm(p_i,iter))])
            disp(['COBE dIdiff: ' num2str(Idiff_max_test_cobe_perm(p_i,iter))])
            disp('-------------------------------------------------------------------------------------------------')
            disp(['Original acco: ' num2str(acco_or_test_perm(p_i,iter))])
            disp(['COBE acco: ' num2str(acco_test_perm(p_i,iter))])

            disp('#################################################################################################')
            save(fullfile(results_path,[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_' res_type '.mat']),"acco_test_perm",'acco_or_test_perm',"Idiff_max_test_perm",...  _cobe_perm_acco_idiff_weight_',num2str(w_idiff*100),'_max2.mat
                "Idiff_max_test_cobe_perm","Idiff_max_test_cobe_idx_perm","cobe_dict",'Component_remv','win','st',"confu_","confu_orr","confu_2","confu_orr2",'Idiff_max_test_cobe_idx_perm_com','Idiff_max_test_cobe_perm_com','-v7.3')
            %     if  acco_test_perm(p_i) > 0.2
            %         disp(acco_test_perm(p_i))
            %     end
        end
    end

    %% Compute the t-statistic
    diff_ = acco_or_test_perm - acco_test_perm;
end

function [acco,acco_or,confu_,confu_orr] = classification_rate_dfc(idf_idx,idf_idx_orr,p,s)

%% Input
% sub_spec --> subject specific FC (train) --> n x p (n = FC upper triangular vector, p = number of training scans)
% sub_spec_test --> subject specific FC (test) --> n x p (p = number of testing scans )
% maph_n --> normalized original FC matrix used to train the dictionary --> n x p1
% maph1_n --> normalized original FC matrix used for testing --> n x p2
% s = sessions per subjects
gt = repmat(1:p,s,1);   % Ground truth
gt = reshape(gt,[],1);   % Ground truth

confu_ = confusionmat(gt,gt(idf_idx));
confu_orr = confusionmat(gt,gt(idf_idx_orr));

acco = sum(diag(confu_))/sum(confu_(:));
acco_or = sum(diag(confu_orr))/sum(confu_orr(:));
end

function [Comm_test,Sub_spec_test] = cobe_transform(vec,COBE_CM_Component)

flat2d_run = reshape(vec,size(vec,1),size(vec,2)*size(vec,3));
X_run_test = COBE_CM_Component'*flat2d_run;
flat2d_Com_run_test = COBE_CM_Component*X_run_test;
flat2d_Sub_spec_run_test = flat2d_run - COBE_CM_Component*X_run_test;

Sub_spec_test = reshape(flat2d_Sub_spec_run_test,size(vec));
Comm_test = reshape(flat2d_Com_run_test,size(vec));
end