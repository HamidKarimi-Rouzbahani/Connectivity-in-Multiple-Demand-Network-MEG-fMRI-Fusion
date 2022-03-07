clc;
clear all;
close all;

subject=1;
Decod_instead_of_Corr=1; % Decoding rates or correlations in RSA? (only for RSA)
Main_analysis_directory='/group/woolgar-lab/projects/Hamid/Projects/MD/Analyses/';
num_iterations=100;
categories=4;
searchlighting=0; % using a searclight(1) or RoI(0)

Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
Beta_folder_name=[Beta_folder_name,'_StmCoding'];



addpath([Main_analysis_directory,'/decoding_toolbox_v3.997/']);
addpath([Main_analysis_directory,'/decoding_toolbox_v3.997/design']);
addpath([Main_analysis_directory,'/decoding_toolbox_v3.997/decoding_results']);
Regions={'ACC_SMA','left_dorsal_lateral_PFc','left_IPS','left_ventral_lateral_PFC',...
    'right_dorsal_lateral_PFc','right_IPS','right_ventral_dorsal_PFc'...
    'LOC_L_10mm_sph_Russ','LOC_R_10mm_sph_Russ','BA1718_2_L_hem','BA1718_1_R_hem'};
addpath(fullfile([Main_analysis_directory,'/spm12']));


csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
csub = csub(1:4);
cfg = decoding_defaults;
cfg.analysis = 'roi';
cfg.decoding.software = 'similarity';
cfg.decoding.method = 'classification';
cfg.decoding.train.classification.model_parameters = 'correlation';
cfg.results.output = 'other';
cfg.plot_selected_voxels = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...
Subject_folder=fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name);
%% Nothing needs to be changed below for a standard similarity analysis using all data
clearvars blocks_with_betas
% finding the blocks containing betas
folders=dir(Subject_folder);
c=0;
for i=1:size(folders,1)
    if strncmp(folders(i).name,'Block',5)
        block_folder_contains=dir(fullfile(Subject_folder,folders(i).name));
        for j=1:size(block_folder_contains,1)
            if strncmp(block_folder_contains(j).name,'beta',4)
                c=c+1;
                blocks_with_betas{c}=folders(i).name;
                break;
            end
        end
    end
end

tmp_file_names=[];
tmp_file_label=[];
tmp_file_descr=[];
tmp_file_chunk=[];
for blk_counter=1:length(blocks_with_betas)
    
    beta_loc = fullfile(Subject_folder,blocks_with_betas{blk_counter},'/');    
    regressor_names = design_from_spm(beta_loc);
    Num_blocks=0;
    for i=1:categories
        Num_blocks=Num_blocks+1;
        labelnames{1,Num_blocks} = regressor_names{1,i};
    end
    labels=[1:categories];
    cfg.results.overwrite = 1;
    cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_loc);
    for i=1:length(cfg.files.name)
        tmp_file_names{(i-1)*length(blocks_with_betas)+blk_counter,1}=cfg.files.name{i};
        tmp_file_label((i-1)*length(blocks_with_betas)+blk_counter,1)=cfg.files.label(i);
        tmp_file_descr{1,(i-1)*length(blocks_with_betas)+blk_counter}=cfg.files.descr{i};
    end
    tmp_file_chunk=vertcat(tmp_file_chunk,blk_counter);
end
cfg.files.name=tmp_file_names;
cfg.files.label=tmp_file_label;
cfg.files.descr=tmp_file_descr;
cfg.files.chunk=repmat(tmp_file_chunk,[max(unique(tmp_file_label)) 1]);


if Decod_instead_of_Corr==0
    cfg.design = make_design_similarity(cfg);
    cfg.basic_checks.DoubleFilenameEntriesOk = 1;
    
    summerization=1; % collapsing across blocks to have one regressor per condition
    if summerization==1
        for class=1:categories
            tmp_beta_cat=0;
            tmp_beta_cat2=0;
            g=0;
            for Num_blocks=[(class-1)*(length(cfg.files.name)./categories)+1:class*(length(cfg.files.name)./categories)]
                g=g+1;
                tmp_beta_cat=tmp_beta_cat+niftiread(cfg.files.name{Num_blocks,1}(1:end));
            end
            
            tmp_beta_cat=tmp_beta_cat./g;
            new_mean_file=fullfile([cfg.files.name{Num_blocks,1}(1:end-13),'mean_c',sprintf('%0.3d', class),'_beta.nii']);
            niftiwrite(tmp_beta_cat,new_mean_file);
            cfg_tmp.files.name{class,1}=new_mean_file;
            cfg_tmp.files.chunk(class,1)=class;
            cfg_tmp.files.label(class,1)=class;
            cfg_tmp.files.descr{1,class}=cfg.files.descr{1,class};
            
        end
        cfg_tmp.files.set=[];
        cfg_tmp.files.xclass=[];
        cfg.files=cfg_tmp.files;
        cfg.design.label=[1:categories]';
        cfg.design.train(categories+1:end)=[];
        cfg.design.test(categories+1:end)=[];
        
    end
    
    if searchlighting==1
        cfg.analysis = 'searchlight';
        cfg.searchlight.unit = 'voxels';
        cfg.searchlight.radious = 100;
        cfg.searchlight.spherical = 1;
        cfg.plot_selected_voxels = 10;
        results = decoding(cfg);
        RSA_results=results;
    else % ROI
        for region=1:length(Regions)
            cfg_tmp_final=cfg;
            for iteration=1:num_iterations
                if iteration>1
                    rand_indd=randsample(1:length(cfg_tmp_final.files.name),length(cfg_tmp_final.files.name));
                    for counter=1:length(cfg_tmp_final.files.name)
                        cfg_tmp_final.files.name{rand_indd(counter),1}=cfg.files.name{counter,1};
                    end
                end
                cfg.files.mask = {[fullfile(Main_analysis_directory, 'Results_temp_playing'),'/ROIs/',csub,'/',Regions{region},'_roi.img']};
                RSA_results{region}=decoding(cfg);
                close all;
            end
        end
    end
    
elseif Decod_instead_of_Corr==1
    
    cfg_tmp=cfg;

    cfg.results.output = {'accuracy','accuracy_minus_chance','balanced_accuracy','sensitivity','specificity','AUC','dprime','loglikelihood','corr','zcorr'}; %, 'SVM_pattern'};
    cfg.decoding.software = 'libsvm';
    cfg.decoding.method = 'classification';
    cfg.analysis = 'roi';
    cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification % parameters for libsvm (linear SV classification, cost = 1, no screen output)
    
    uniqe_labels=unique(cfg.files.label);
    combinations=nchoosek([1:length(uniqe_labels)],2);
    for comb=1:length(combinations)
        ind1=(cfg_tmp.files.label==combinations(comb,1));
            ind2=(cfg_tmp.files.label==combinations(comb,2));

        t=0;
        for n=1:length(ind1)
            if ind1(n)==1
                t=t+1;
                name_tmp{t,1}=cfg_tmp.files.name{n};
                chunk_tmp(t,1)=cfg_tmp.files.chunk(n);
                label_tmp(t,1)=cfg_tmp.files.label(n);
                desc_tmp{1,t}=cfg_tmp.files.descr{1,n};
            end
        end
        for n=1:length(ind2)
            if ind2(n)==1
                t=t+1;
                name_tmp{t,1}=cfg_tmp.files.name{n};
                chunk_tmp(t,1)=cfg_tmp.files.chunk(n);
                label_tmp(t,1)=cfg_tmp.files.label(n);
                desc_tmp{1,t}=cfg_tmp.files.descr{1,n};
            end
        end

        
        cfg.files.name=name_tmp;
        cfg.files.chunk=chunk_tmp;
        cfg.files.label=label_tmp;
        cfg.files.descr=desc_tmp;
        %% equalization of samples in training
        unq_labls=unique(cfg.files.label);
        for i=1:length(unq_labls)
            inds_lab(i)=sum(cfg.files.label==unq_labls(i));
        end
        if abs(diff(inds_lab))>0
            [min_num,which_class]=min(inds_lab);
            [inds_cls_to_remove_from,~]=find(cfg.files.label==unq_labls(length(unq_labls)+1-which_class));
            ind_to_remove=randsample(inds_cls_to_remove_from,abs(diff(inds_lab)));
            cfg.files.name(ind_to_remove)=[];
            cfg.files.chunk=[[1:min(inds_lab)]';[1:min(inds_lab)]'];
            cfg.files.label(ind_to_remove)=[];
            cfg.files.descr(ind_to_remove)=[];
        end
        %% CV and decoding
        cfg.design = make_design_cv(cfg);
        cfg.basic_checks.DoubleFilenameEntriesOk = 1;
        
        if searchlighting==1
            cfg.analysis = 'searchlight';
            cfg.searchlight.unit = 'voxels';
            cfg.searchlight.radious = 100;
            cfg.searchlight.spherical = 1;
            cfg.plot_selected_voxels = 100;
            results = decoding(cfg);
            RSA_results{comb}=results;
        else
            for region=1:length(Regions)
                cfg_tmp_final=cfg;
                for iteration=1:num_iterations
                    if iteration>1
                        rand_indd=randsample(1:length(cfg_tmp_final.files.name),length(cfg_tmp_final.files.name));
                        for counter=1:length(cfg_tmp_final.files.name)
                            cfg_tmp_final.files.name{rand_indd(counter),1}=cfg.files.name{counter,1};
                        end
                    end
                    cfg_tmp_final.files.mask = {[fullfile(Main_analysis_directory, 'Results_temp_playing'),'/ROIs/',csub,'/',Regions{region},'_roi.img']};
                    RSA_results{region,comb,iteration}=decoding(cfg_tmp_final);
                    close all;
                end
            end
        end
    end
end
Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];
if Decod_instead_of_Corr==1
    dec='Dec_';
else
    dec='';
end
if num_iterations>1
    Iter='Iter_';
else
    Iter='';
end

save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'CBU_Stim_side_acrs_rules',dec,SL,'.mat']),'RSA_results');

