% close all
clear
working_dir = 'H:\Single_scan_cobe_results2';

addpath(fullfile(working_dir,'codes_functions','functions'))
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks


atlas_name = 's400_hcp'; %s400_hcp a116

% full atlas name
switch atlas_name
    case 's400_hcp'
        atlas_full_name = 'Schaefer-HCP 419';
        net_select_perm = 1:9; 
    case 'a116'
        atlas_full_name = 'AAL 116';
        net_select_perm = [1:3,5:9];
    case 'w300'
        atlas_full_name = 'Seitzman 300';
        net_select_perm = 1:9;
    case 'p264'
        atlas_full_name = 'Power 264';
        net_select_perm = 1:9;
    case 'd164'
        atlas_full_name = 'Dosenbach 164';
        net_select_perm = [1:5,7:9];
end


% net_select = 1:8;

Win_size = [50,100,200,400,800];
Com_comp = 1:2:20;
Stride = [20,50,100];
Stride_s = [14.4, 36, 72];
stride_size = length(Stride);


%% Parameter check
results_path = fullfile(working_dir,'para_check');
net_select = 9;%1:net_size;
figure;
for sub_i = 1:3

    for i_net = net_select
        load(fullfile(results_path,[atlas_name,'_para_test',num2str(sub_i),'_',char(net_names{i_net}),'.mat'])) %_with_new_windows
        if size(acco,2) == 4
            acco = acco(:,2:end,:);
            acco_or = acco_or(:,2:end,:);
        end
        acco(acco==0) = nan;
        acco = permute(acco,[1,3,2]);
        [~,com_size,stride_size] = size(acco);

        for i_s = 1:stride_size
            subplot(3,3,(sub_i-1)*3 + i_s)
            plot((acco(:,:,i_s)'*100),'LineWidth',2)
            title(['dFC Stride = ' num2str(Stride_s(i_s))  ' s' ', Block = ' num2str(sub_i)])
            xlabel('Number of Common COBE components')
            ylabel('Identification Rate')
            %         set(gca,'ytick',1:win_size,'yticklabel',Win)
            set(gca,'xtick',1:com_size,'xticklabel',Com_comp)
            %         yticks(Win)
            %         xticks(Com_comp)
            grid on
            set(gca,'FontSize',16)
            ylim([60 110])
        end
    end
    
end
legend({'36 s','72 s','144 s','288 s','576 s'}) %'36 s','72 s','144 s','288 s',  ,'720 s','792 s','828 s'
%%
%% Permutations check
% di_diff    
results_path = fullfile(working_dir,'5x2cv',atlas_name);
figure;

for i_net = net_select_perm 
    load(fullfile(results_path,[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2.mat']),'Idiff_max_test_perm','Idiff_max_test_cobe_perm','idiff_score_s')
%     net_acco_test_s(:,i_net) = idiff_score_s(:);
    net_acco_test(:,i_net) = Idiff_max_test_perm(:)*100;

    net_acco_test_cobe(:,i_net) = Idiff_max_test_cobe_perm(:)*100;

%     [~,p_idiff(i_net)] = ttest2(net_acco_test_cobe(:,i_net),net_acco_test(:,i_net));

    [t(i_net),p(i_net)] = p5x2cv(Idiff_max_test_cobe_perm,Idiff_max_test_perm);
end
p= p';
h = p<0.0056;
net_all = cat(3,net_acco_test,net_acco_test_cobe);%cat(3,net_acco_test_s,net_acco_test,net_acco_test_cobe);


barplot_with_errorbars(net_all,net_names);ylim([0 80])
legend({'Original dFC','COBE Subject-Specific dFC'})
set(gcf,"Position",[1,1,1850,500])
ylabel('dI_{diff} values')
xlabel('Resting state Networks')
set(gca,'FontSize',17)
title(atlas_full_name)
grid on

%% acco
  
figure;
clear net_acco_test  net_acco_test_cobe
for i_net = net_select_perm 
    load(fullfile(results_path,[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2.mat']),'acco_or_test_perm','acco_test_perm')
%     net_acco_test_s(:,i_net) = acco_or_test_perm_s(:)*100;
    net_acco_test(:,i_net) = acco_or_test_perm(:)*100;

    net_acco_test_cobe(:,i_net) = acco_test_perm(:)*100;

%     [~,p_acco(i_net)] = ttest2(net_acco_test_cobe(:,i_net),net_acco_test(:,i_net));
    [t(i_net),p(i_net)] = p5x2cv(acco_test_perm,acco_or_test_perm);
end
net_all = cat(3,net_acco_test,net_acco_test_cobe); % net_all = cat(3,net_acco_test_s,net_acco_test,net_acco_test_cobe); %
p= p';
h = p<0.0056;

barplot_with_errorbars(net_all,net_names);ylim([0 120])
legend({'Original dFC','COBE-Subject-Specific dFC'}) %legend({'sFC','Original dFC','COBE-Subject-Specific dFC'})
title([])
ylabel('Identification Rate')
xlabel('Resting state Networks')
set(gcf,"Position",[1,1,1850,500])
set(gca,'FontSize',17)
title(atlas_full_name)
grid on
%% idiff with family data

 
results_path = fullfile(working_dir,'5x2cv_twins',atlas_name);
figure;
clear net_acco_test  net_acco_test_cobe
for i_net = net_select_perm 
    load(fullfile(results_path,'nt_un',[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_nt_un.mat']),'Idiff_max_test_perm','Idiff_max_test_cobe_perm','idiff_score_s')

    net_acco_test_un(:,i_net) = Idiff_max_test_perm(:)*100;

    net_acco_test_cobe_un(:,i_net) = Idiff_max_test_cobe_perm(:)*100;
    
    load(fullfile(results_path,'MZ',[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_MZ.mat']),'Idiff_max_test_perm','Idiff_max_test_cobe_perm','idiff_score_s')
    net_acco_test_mz(:,i_net) = Idiff_max_test_perm(:)*100;

    net_acco_test_cobe_mz(:,i_net) = Idiff_max_test_cobe_perm(:)*100;
    
    load(fullfile(results_path,'DZ',[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_DZ.mat']),'Idiff_max_test_perm','Idiff_max_test_cobe_perm','idiff_score_s')
    net_acco_test_dz(:,i_net) = Idiff_max_test_perm(:)*100;

    net_acco_test_cobe_dz(:,i_net) = Idiff_max_test_cobe_perm(:)*100;

end
net_all = cat(3,net_acco_test_mz,net_acco_test_dz,net_acco_test_un,net_acco_test_cobe_mz,net_acco_test_cobe_dz,net_acco_test_cobe_un); % net_all = cat(3,net_acco_test_s,net_acco_test,net_acco_test_cobe); %
barplot_with_errorbars(net_all,net_names);ylim([0 120])
legend({'Original dFC MZ','Original dFC DZ','Original dFC un','COBE-Subject-Specific dFC MZ','COBE-Subject-Specific dFC DZ','COBE-Subject-Specific dFC un'}) %legend({'sFC','Original dFC','COBE-Subject-Specific dFC'})
% legend({'Original sFC','Original dFC','Original sFC MZ','Original dFC MZ','Original sFC un','Original dFC un','COBE-Subject-Specific dFC','COBE-Subject-Specific dFC MZ','COBE-Subject-Specific dFC un'}) %legend({'sFC','Original dFC','COBE-Subject-Specific dFC'})

title([])
ylabel('dI_{diff}')
xlabel('Resting state Networks')
set(gcf,"Position",[1,1,1850,500])
set(gca,'FontSize',17)
title(atlas_full_name)
grid on

%% acco with family data
%  acco
figure;
clear net_acco_test  net_acco_test_cobe
for i_net = net_select_perm 
    load(fullfile(results_path,'nt_un',[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_nt_un.mat']),'acco_or_test_perm_s','acco_or_test_perm','acco_test_perm')

    net_acco_test_un(:,i_net) = acco_or_test_perm(:)*100;

    net_acco_test_cobe_un(:,i_net) = acco_test_perm(:)*100;
    
    load(fullfile(results_path,'MZ',[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_MZ.mat']),'acco_or_test_perm_s','acco_or_test_perm','acco_test_perm')
    net_acco_test_mz(:,i_net) = acco_or_test_perm(:)*100;

    net_acco_test_cobe_mz(:,i_net) = acco_test_perm(:)*100;
    
    load(fullfile(results_path,'DZ',[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2_DZ.mat']),'acco_or_test_perm_s','acco_or_test_perm','acco_test_perm')
    net_acco_test_dz(:,i_net) = acco_or_test_perm(:)*100;

    net_acco_test_cobe_dz(:,i_net) = acco_test_perm(:)*100;
    
    [t_mz_un(i_net),p_mz_un(i_net)] = p5x2cv(reshape(net_acco_test_mz(:,i_net),2,5),reshape(net_acco_test_un(:,i_net),2,5));
    [t_dz_un(i_net),p_dz_un(i_net)] = p5x2cv(reshape(net_acco_test_dz(:,i_net),2,5),reshape(net_acco_test_un(:,i_net),2,5));
    [t_mz_dz(i_net),p_mz_dz(i_net)] = p5x2cv(reshape(net_acco_test_mz(:,i_net),2,5),reshape(net_acco_test_dz(:,i_net),2,5));

    [t_mz_c(i_net),p_mz_c(i_net)] = p5x2cv(reshape(net_acco_test_mz(:,i_net),2,5),reshape(net_acco_test_cobe_mz(:,i_net),2,5));
    [t_un_c(i_net),p_un_c(i_net)] = p5x2cv(reshape(net_acco_test_un(:,i_net),2,5),reshape(net_acco_test_cobe_un(:,i_net),2,5));

end
net_all = cat(3,net_acco_test_mz,net_acco_test_dz,net_acco_test_un,net_acco_test_cobe_mz,net_acco_test_cobe_dz,net_acco_test_cobe_un); % net_all = cat(3,net_acco_test_s,net_acco_test,net_acco_test_cobe); %
barplot_with_errorbars(net_all,net_names);ylim([0 120])
legend({'Original dFC MZ','Original dFC DZ','Original dFC un','COBE-Subject-Specific dFC MZ','COBE-Subject-Specific dFC DZ','COBE-Subject-Specific dFC un'}) %legend({'sFC','Original dFC','COBE-Subject-Specific dFC'})

title([])
ylabel('Identification Rate')
xlabel('Resting state Networks')
set(gcf,"Position",[1,1,1850,500])
set(gca,'FontSize',17)
title(atlas_full_name)
grid on

%% COBE Dictionary (more than 1 basis)
results_path = fullfile(working_dir,'5x2cv',atlas_name);
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks

i_net = 9;
atlas_path = fullfile(working_dir,'atlases',[atlas_name,'.nii']);
atlas_yeo = fullfile(working_dir,'atlases','yeo_7_net.nii');
idx = yeo_networks(atlas_path,atlas_yeo);
[sort_idx_val,sort_idx] = sort(idx);
load(fullfile(results_path,[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2.mat'])) %[char(net_names(i_net)),'_max_acco1_' atlas_name '_test5x2.mat']
% figure;imagesc(reshape(cobe_dict,size(cobe_dict,1),[]))
figure;imagesc(abs(corr(reshape(cobe_dict,size(cobe_dict,1),[]))));axis image;title('Correlation of all Dictionary basis')
cobe_dic_map =veccorr(squeeze(cobe_dict(:,1,1)));

% Average value per net connection
dict_map = cobe_dic_map(sort_idx,sort_idx);
sort_idx_net = sort(idx);
%
net_sum_cobe_dict = zeros(8);
net_avg_cobe_dict = zeros(8);
for i = 1:8
    for j = 1:8
       net_sum_cobe_dict(i,j) = sum(dict_map(sort_idx_net==i-1,sort_idx_net==j-1),'all');
       net_avg_cobe_dict(i,j) = net_sum_cobe_dict(i,j)/(sum(nnz(sort_idx_net==i-1))*sum(nnz(sort_idx_net==j-1)));
    end
end

figure;
subplot(1,2,1)
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode'}; % names of the networks
fc_imshow(cobe_dic_map(sort_idx,sort_idx),sort_idx_val,net_names);axis image  ;%clim([-0.015 0.01])%title('COBE Most reproducible basis')
set(gca,'FontSize',15)
clim([-10 5]*10^(-3))
% colormap("jet")
subplot(1,2,2)
imagesc(net_avg_cobe_dict);axis image;%clim([-6 3]*10^(-3))
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode'}; % names of the networks
set(gca,'xtick',1:8,'xticklabel',net_names)
set(gca,'ytick',1:8,'yticklabel',net_names)
set(gca,'FontSize',15)
clim([-0.01 0])
colorbar
% colormap("jet")

%% From HCP dict, NKI test


%% Permutations check
% di_diff    
net_names = {'Non Yeo','Visual','Somato Motor','Dorsal Attention','Ventral Attention','Limbic Network','Fronto-Parietal','Default Mode','Whole Brain'}; % names of the networks

results_path = ['C:\Users\krish\Box\Single Subject COBE\new_data_results\nki_hcp\from_hcp_dict\new_acco\' atlas_name];
figure;
clear net_acco_test_cobe net_acco_test
for i_net = net_select_perm 
    load(fullfile(results_path,[char(net_names(i_net)),'_nki_' atlas_name '_test_from_hcp_perm_dict_new_acco_5x2.mat']),'Idiff_max_test_perm_nki','Idiff_max_test_cobe_perm_nki')
    net_acco_test(:,i_net) = Idiff_max_test_perm_nki(:)*100;

    net_acco_test_cobe(:,i_net) = Idiff_max_test_cobe_perm_nki(:)*100;

    [t(i_net),p(i_net)] = p5x2cv(net_acco_test_cobe(:,i_net),net_acco_test(:,i_net));
end
p = p';
h = p<0.0056;
net_all = cat(3,net_acco_test,net_acco_test_cobe);% cat(3,net_acco_test_s,net_acco_test,net_acco_test_cobe);


barplot_with_errorbars(net_all,net_names);ylim([0 80])
legend({'Original dFC','COBE Subject-Specific dFC'})
set(gcf,"Position",[1,1,1850,600])
ylabel('dI_{diff} values')
xlabel('Resting state Networks')
set(gca,'FontSize',17)
title([atlas_full_name ', NKI Subject-Specific, COBE dictionaries trained using HCP'])
grid on

%% acco
results_path = fullfile(working_dir,'from_hcp_dict_nki',atlas_name);
figure;
clear net_acco_test  net_acco_test_cobe
for i_net = net_select_perm 
    load(fullfile(results_path,[char(net_names(i_net)),'_nki_' atlas_name '_test_from_hcp_perm_dict_5x2.mat']),'acco_or_test_perm_nki','acco_test_perm_nki')
%     net_acco_test_s(:,i_net) = acco_or_test_perm_s(:)*100;
    net_acco_test(:,i_net) = acco_or_test_perm_nki(:)*100;

    net_acco_test_cobe(:,i_net) = acco_test_perm_nki(:)*100;

     [t(i_net),p(i_net)] = p5x2cv(net_acco_test_cobe(:,i_net),net_acco_test(:,i_net));
end
p = p';
h = p<0.0056;
net_all = cat(3,net_acco_test,net_acco_test_cobe); %net_all = cat(3,net_acco_test_s,net_acco_test,net_acco_test_cobe);


barplot_with_errorbars(net_all,net_names);ylim([0 120])
legend({'Original dFC','COBE-Subject-Specific dFC'})
title([])
ylabel('Identification Rate')
xlabel('Resting state Networks')
set(gcf,"Position",[1,1,1850,600])
set(gca,'FontSize',17)
title([atlas_full_name ', NKI Subject-Specific, COBE dictionaries trained using HCP'])
grid on


