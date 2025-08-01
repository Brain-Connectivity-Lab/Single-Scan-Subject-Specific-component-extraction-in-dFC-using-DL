%% Create 3 blocks and reg_ts
% This script was written by Pratik Jain if you have any questions email
% jainpratik412@gmail.com
% This Script creates the average ROI timeseries for all the subjects and
% divides the subjects into 3 blocks randomly.
% We ran this for both the LR and RL encodings
% then 

%% Data and Results path
data_folder = 'H:\HCP_S1200_preprocessed_denoised';
working_dir = 'H:\Single_scan_cobe_results2';

%% Atlas
atlas_path = fullfile(working_dir,'atlases');
at_names = {'a116','d164','p264','w300','s400_hcp'}; %AAL, Dosenbach, Power, Seitzman, Shaefer_HCP 

%% add paths

addpath(fullfile(working_dir,'codes_functions','functions'))

run = '1';
encode= 'LR';
all_sub_files_LR = dir(fullfile(data_folder,'*',['MNINonLinear\Results\rfMRI_REST' run '_' encode],['*_s400_run_' run '_encode_' encode '_reg_ts.mat']));
encode= 'RL';
all_sub_files_RL = dir(fullfile(data_folder,'*',['MNINonLinear\Results\rfMRI_REST' run '_' encode],['*_s400_run_' run '_encode_' encode '_reg_ts.mat']));

LR_names = struct2cell(all_sub_files_LR); LR_names(2:end,:) = []; LR_names = cellfun(@(x) x(1:6),LR_names,'UniformOutput',false);
RL_names = struct2cell(all_sub_files_RL); RL_names(2:end,:) = []; RL_names = cellfun(@(x) x(1:6),RL_names,'UniformOutput',false);

[~,ai,bi] = intersect(string(LR_names),string(RL_names));

all_sub_files_LR = all_sub_files_LR(ai);
all_sub_files_RL = all_sub_files_RL(bi);

matfile_path = fullfile(working_dir,'matfiles'); 
mkdir(matfile_path)

n_subj = length(all_sub_files_LR);

n_sub_block = round(n_subj/3);

all_sub = struct2cell(all_sub_files_RL);
all_sub(2:end,:) = [];
for i = 1:length(all_sub)
    all_sub{i} = all_sub{i}(1:6);
end

subj_idx_randperm = randperm(n_subj);

sub_block_idx{1} = subj_idx_randperm(1:n_sub_block);
sub_block_idx{2} = subj_idx_randperm(n_sub_block+1:2*n_sub_block);
sub_block_idx{3} = subj_idx_randperm(2*n_sub_block+1:end);

for at_i = 5     % Atlas for loop put the Atlas number that you want to create the matfile for here
    atlas = niftiread(fullfile(atlas_path,[char(at_names{at_i}),'.nii']));
    level = unique(sort((atlas(:))));
    level = level(2:end);

    for i = 1:length(sub_block_idx)
        p_sub_name = cell(length(sub_block_idx{i}),1);
        reg_ts = zeros(1200,length(level(2:end)),length(sub_block_idx{i}));
        vox_count_in_reg = zeros(length(level(2:end)),length(sub_block_idx{i}));
        for i_sub = 1:length(sub_block_idx{i})

            p_sub_name{i_sub} = all_sub_files_RL(sub_block_idx{i}(i_sub)).name(1:6);

            if strcmp(encode,'LR')
                subject = niftiread(complete_filepath(all_sub_files_LR(sub_block_idx{i}(i_sub)).folder,'rfMRI*_hp2000_clean.nii.gz'));
            else
                subject = niftiread(complete_filepath(all_sub_files_RL(sub_block_idx{i}(i_sub)).folder,'rfMRI*_hp2000_clean.nii.gz'));
            end


            [x,y,z,nT] = size(subject);
            subMat = reshape(subject,[x*y*z,nT]);
            subMat = double(subMat');
            subMat_mean = mean(subMat,1);


            for r = 1:length(level)
                vox_id = find(atlas==level(r));
                reg_voxels_ts = subMat(:,vox_id) - subMat_mean(vox_id);
                vox_count_in_reg(r,i_sub) = size(reg_voxels_ts,2);
                reg_ts(:,r,i_sub) = mean(reg_voxels_ts,2);  % mean TS of region 'r'
            end

            disp(['Subject ' p_sub_name{i_sub} ' from sub-block ' num2str(i) ' done. Number ' num2str(i_sub) ])
            save(fullfile(matfile_path,['HCP_',char(at_names{at_i}),'_run_',run,'_encode_',encode,'_sub_block',num2str(i),'_reg_ts.mat']),'reg_ts','p_sub_name','vox_count_in_reg')
        end
        save(fullfile(matfile_path,['HCP_',char(at_names{at_i}),'_run_',run,'_encode_',encode,'_sub_block',num2str(i),'_reg_ts.mat']),'reg_ts','p_sub_name','vox_count_in_reg')
    end

end