function [done]=Step8_Decoding_f_pBlk3(subject,Main_analysis_directory,Action,Trial_type,Decod_instead_of_Corr,searchlighting)

num_iterations=101;
if Action==1
    Beta_folder_name='Correct_trials_decoding_design_pBlk_native';
elseif Action==2
    Beta_folder_name='Correct_trials_RSA_design_16_pBlk_native';
    categories=16;
elseif Action==3
    Beta_folder_name='Correct_trials_RSA_design_8_pBlk_native';
    categories=8;
elseif Action==4
    Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
    categories=4;
end
if Trial_type==2
    Beta_folder_name_prim=[Beta_folder_name];
    Beta_folder_name_scnd=[Beta_folder_name,'_AllErr'];
elseif Trial_type==3
    if Action==1
        Beta_folder_name_prim=[Beta_folder_name];
        Beta_folder_name_scnd=[Beta_folder_name,'_OpsErr'];
    elseif Action>1
        Beta_folder_name=[Beta_folder_name,'_CueCoding'];
    end
elseif Trial_type==4
    if Action==1
        Beta_folder_name=[Beta_folder_name,'_IncldSide'];
    elseif Action>1
        Beta_folder_name=[Beta_folder_name,'_StmCoding'];
    end
elseif Trial_type==5
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_StmPeriCod_clpsd_OppstCuCol'];
    end
elseif Trial_type==6
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_StmSideCod_clpsd_OppstCuCol'];
    end
elseif Trial_type==7
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_RespCod_1to4'];
    end
elseif Trial_type==8
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_RespCod_inRules'];
    end
elseif Trial_type==9
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_RespCod_clpsd_OppstCuCol'];
    end
elseif Trial_type==10
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_StimSideCoding'];
    end
elseif Trial_type==11
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_RespCod_InOut_clpsd_OppstCuCol'];
    end
elseif Trial_type==12
    if Action>1
        Beta_folder_name=[Beta_folder_name,'_RespCod_InOut_clpsd2_OppstCuCol'];
    end
end



addpath([Main_analysis_directory,'/decoding_toolbox_v3.997/']);
addpath([Main_analysis_directory,'/decoding_toolbox_v3.997/design']);
addpath([Main_analysis_directory,'/decoding_toolbox_v3.997/decoding_results']);
Regions={'ACC_SMA','left_dorsal_lateral_PFc','left_IPS','left_ventral_lateral_PFC',...
    'right_dorsal_lateral_PFc','right_IPS','right_ventral_dorsal_PFc'...
    'LOC_L_10mm_sph_Russ','LOC_R_10mm_sph_Russ','BA1718_2_L_hem','BA1718_1_R_hem'};
addpath(fullfile([Main_analysis_directory,'\spm12']));

if Action==1
    if Trial_type==1
        decodings=1:2;
    elseif Trial_type==4
        decodings=1:4;
    end
    for decod=decodings
        
        clearvars -except Decod_instead_of_Corr searchlighting decod Beta_folder_name_scnd Beta_folder_name_prim categories Regions Beta_folder_name num_iterations subject Action Trial_type Behavioural_directory Behavioural_directory fmri_data_address Main_analysis_directory
        
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        
        cfg = decoding_defaults;
        cfg.results.output = {'accuracy','accuracy_minus_chance','balanced_accuracy','sensitivity','specificity','AUC','dprime','loglikelihood','corr','zcorr'}; %, 'SVM_pattern'};
        
        cfg.decoding.method = 'classification';
        cfg.analysis = 'roi';
        
        cfg.results.write = 1;
        cfg.results.overwrite = 1;
        cfg.plot_selected_voxels = 0;
        cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification % parameters for libsvm (linear SV classification, cost = 1, no screen output)
        cfg.scale.method = 'none';
        cfg.scale.estimation = 'none';
        
        Subject_folder=fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name);
        
        
        Num_blocks=0;
        if decod==1
            classes={'S1C','S2C'};
        elseif decod==2
            classes={'R1C','R2C'};
        elseif decod==3
            classes={'P1C','P2C'};
        elseif decod==4
            classes={'L1C','L2C'};
        end
        
        
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
            temp_block_folder=fullfile(Subject_folder,blocks_with_betas{blk_counter},'/');
            regressor_names = design_from_spm(temp_block_folder);
            cfg = decoding_describe_data(cfg,classes,[1 2],regressor_names,temp_block_folder);
            
            tmp_file_names{blk_counter,1}=cfg.files.name{1};
            tmp_file_names{blk_counter+length(blocks_with_betas),1}=cfg.files.name{2};
            tmp_file_label(blk_counter,1)=cfg.files.label(1);
            tmp_file_label(blk_counter+length(blocks_with_betas),1)=cfg.files.label(2);
            tmp_file_descr{1,blk_counter}=cfg.files.descr{1};
            tmp_file_descr{1,blk_counter+length(blocks_with_betas)}=cfg.files.descr{2};
            tmp_file_chunk=vertcat(tmp_file_chunk,blk_counter);
            
        end
        cfg.files.name=tmp_file_names;
        cfg.files.label=tmp_file_label;
        cfg.files.descr=tmp_file_descr;
        cfg.files.chunk=repmat(tmp_file_chunk,[2 1]);
        cfg.design = make_design_cv(cfg);
        
        
        if Trial_type==2 || Trial_type==3
            if Trial_type==2
                if decod==1
                    classes={'S1S','S2S'};
                elseif decod==2
                    classes={'R1R','R2R'};
                elseif decod==3
                    classes={'P1P','P2P'};
                end
            elseif Trial_type==3
                if decod==1
                    classes={'S1R','S2R'};
                elseif decod==2
                    classes={'R1S','R2S'};
                elseif decod==3
                    classes={'P1E','P2E'};
                end
            end
            
            beta_loc_secondary = fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name_scnd);
            clearvars blocks_with_betas
            % finding the blocks containing betas
            folders=dir(beta_loc_secondary);
            c=0;
            for i=1:size(folders,1)
                if strncmp(folders(i).name,'Block',5)
                    block_folder_contains=dir(fullfile(beta_loc_secondary,folders(i).name));
                    for j=1:size(block_folder_contains,1)
                        if strncmp(block_folder_contains(j).name,'beta',4)
                            c=c+1;
                            blocks_with_betas{c}=folders(i).name;
                            break;
                        end
                    end
                end
            end
            
            cfg.design.train=repmat(cfg.design.train,[1 size(cfg.design.train,2)]);
            cfg.design.train=vertcat(cfg.design.train,zeros(length(blocks_with_betas)*2,size(cfg.design.train,2)));
            
            tmp_file_names=[];
            tmp_file_label=[];
            tmp_file_descr=[];
            tmp_file_chunk=[];
            for blk_counter=1:length(blocks_with_betas)
                temp_block_folder=fullfile(beta_loc_secondary,blocks_with_betas{blk_counter},'/');
                regressor_names = design_from_spm(temp_block_folder);
                cfg2 = decoding_describe_data(cfg,classes,[1 2],regressor_names,temp_block_folder);
                tmp_file_names{blk_counter,1}=cfg2.files.name{1};
                tmp_file_names{blk_counter+length(blocks_with_betas),1}=cfg2.files.name{2};
                tmp_file_label(blk_counter,1)=cfg2.files.label(1);
                tmp_file_label(blk_counter+length(blocks_with_betas),1)=cfg2.files.label(2);
                tmp_file_descr{1,blk_counter}=cfg2.files.descr{1};
                tmp_file_descr{1,blk_counter+length(blocks_with_betas)}=cfg2.files.descr{2};
                tmp_file_chunk=vertcat(tmp_file_chunk,blk_counter);
            end
            cfg2.files.name=tmp_file_names;
            cfg2.files.label=tmp_file_label;
            cfg2.files.descr=tmp_file_descr;
            cfg2.files.chunk=repmat(tmp_file_chunk,[2 1]);
            cfg2.design = make_design_cv(cfg2);
            
            % all error trials tested in every cross-validation round
            cfg2.design.test=zeros(size(cfg.design.train));
            cfg2.design.test(size(cfg.design.test,1)+1:end,:)=ones(size(cfg2.design.train,1),size(cfg2.design.test,2));
            
            cfg.design.test=cfg2.design.test;
            cfg.design.set=repmat(cfg.design.set,[1 length(cfg.design.set)]);
            cfg.files.label=vertcat(cfg.files.label,cfg2.files.label);
            cfg.design.label=repmat(cfg.files.label,[1 length(cfg.design.set)]);
            cfg.files.name=vertcat(cfg.files.name,cfg2.files.name);
            cfg.files.chunk=vertcat(cfg.files.chunk,cfg2.files.chunk);
            cfg.files.descr=horzcat(cfg.files.descr,cfg2.files.descr);
            
            % remove empty (no regressor) beta files and equalize betas across conditions
            for beta_number=1:size(cfg.files.name,1)
                if nanmean(nanmean(nanmean(niftiread(cfg.files.name{beta_number,1}(1:end)))))==0
                    empty_betas(beta_number)=1;
                else
                    empty_betas(beta_number)=0;
                end
            end
            
            for cat=1:2
                num_empty_betas_per_cat(cat)=sum(cfg.files.label(empty_betas==1)==cat);
            end
            [~,cat_with_more_samples]=min(num_empty_betas_per_cat);
            num_samples_to_remove=abs(diff(num_empty_betas_per_cat));
            second_half_samples=([length(cfg.files.label)-length(cfg2.files.label)+1:length(cfg.files.label)]);
            ind_other_class_to_remove=find([cfg.files.label(second_half_samples)==cat_with_more_samples & [empty_betas(second_half_samples)==0]']);
            empty_betas(length(cfg.files.label)-length(cfg2.files.label)+randsample(ind_other_class_to_remove,num_samples_to_remove))=1;
            
            c=0;
            for beta_number=1:size(cfg.files.name,1)
                if empty_betas(beta_number)==0
                    c=c+1;
                    label_remained(c,:)=cfg.design.label(beta_number,:);
                    train_remained(c,:)=cfg.design.train(beta_number,:);
                    test_remained(c,:)=cfg.design.test(beta_number,:);
                    
                    name_remained{c,1}=cfg.files.name{beta_number,1};
                    label_files_remained(c,1)=cfg.files.label(beta_number);
                    descr_remained{1,c}=cfg.files.descr{1,beta_number};
                    chunk_remained(c,1)=cfg.files.chunk(beta_number);
                end
            end
            cfg.design.label=label_remained;
            cfg.design.train=train_remained;
            cfg.design.test=test_remained;
            cfg.files.name=name_remained;
            cfg.files.label=label_files_remained;
            cfg.files.descr=descr_remained;
            cfg.files.chunk=chunk_remained;
        elseif Trial_type==4
            if decod==1
                classes={'S1C','S2C'};
            elseif decod==2
                classes={'R1C','R2C'};
            elseif decod==3
                classes={'P1C','P2C'};
            elseif decod==4
                classes={'L1C','L2C'};
            end
        end
        if searchlighting==1
            cfg.analysis = 'searchlight';
            cfg.searchlight.unit = 'voxels';
            cfg.searchlight.radious = 50;
            cfg.searchlight.spherical = 1;
            cfg.plot_selected_voxels = 10;
            results = decoding(cfg);
            Decoding_results=results;
        else
            cfg_saved=cfg;
            for region=1:length(Regions)
                cfg=cfg_saved;
                cfg.files.mask = {[fullfile(Main_analysis_directory, 'Results_temp_playing'),'\ROIs\',csub,'\',Regions{region},'_roi.img']};
                
                results = decoding(cfg);
                Decoding_results{region,1}=results;
                for iteration=2:num_iterations
                    tmp_labels=nan*ones(size(cfg.design.label,2)*2,size(cfg.design.label,2));
                    for col=1:size(cfg.design.label,2)
                        tmp_labels([col col+size(cfg.design.label,2)],col)=[1;2];
                        other_labels=randsample([ones(size(cfg.design.label,2)-1,1);2*ones(size(cfg.design.label,2)-1,1)],size(cfg.design.label,2)*2-2);
                        [inds_nans,~]=find(isnan(tmp_labels(:,col)));
                        tmp_labels(inds_nans,col)=other_labels;
                        cfg.design.label(:,col)=tmp_labels(:,col);
                    end
                    results = decoding(cfg);
                    Decoding_results{region,iteration}=results;
                    close all;
                end
                close all;
            end
        end
        if searchlighting==1
            SL='SearchLight';
        else
            SL='ROIs';
        end
        
        Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];
        if Trial_type==2 || Trial_type==3
            save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name_scnd,'/',['Dec_',SL,'_',classes{1,1}(1),'.mat']),'Decoding_results');
        elseif Trial_type==1 || Trial_type==4
            save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_',SL,'_',classes{1,1}(1),'.mat']),'Decoding_results');
        end
    end
elseif Action>1
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
        if Trial_type==2
            beta_loc = fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name_prim,blocks_with_betas{blk_counter},'/');
            beta_loc_secondary = fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name_scnd,blocks_with_betas{blk_counter},'/');
        end
        %         cfg.results.dir = beta_loc;
        
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
    
    if Trial_type==2
        % finding the blocks containing betas
        beta_loc_secondary = fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name_scnd);
        folders=dir(beta_loc_secondary);
        clearvars blocks_with_betas
        c=0;
        for i=1:size(folders,1)
            if strncmp(folders(i).name,'Block',5)
                block_folder_contains=dir(fullfile(beta_loc_secondary,folders(i).name));
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
            beta_loc_secondary = fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name_scnd,blocks_with_betas{blk_counter},'/');
            
            %                 cfg.results.dir = beta_loc;
            
            regressor_names = design_from_spm(beta_loc);
            Num_blocks=0;
            for i=1:categories
                Num_blocks=Num_blocks+1;
                labelnames{1,Num_blocks} = regressor_names{1,i};
            end
            labels=[1:categories];
            cfg.results.overwrite = 1;
            cfg2 = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_loc_secondary);
            
            for i=1:length(cfg2.files.name)
                tmp_file_names{(i-1)*length(blocks_with_betas)+blk_counter,1}=cfg2.files.name{i};
                tmp_file_label((i-1)*length(blocks_with_betas)+blk_counter,1)=cfg2.files.label(i);
                tmp_file_descr{1,(i-1)*length(blocks_with_betas)+blk_counter}=cfg2.files.descr{i};
            end
            tmp_file_chunk=vertcat(tmp_file_chunk,blk_counter);
        end
        cfg2.files.name=tmp_file_names;
        cfg2.files.label=tmp_file_label;
        cfg2.files.descr=tmp_file_descr;
        %             cfg2.files.chunk=tmp_file_chunk;
        cfg2.files.chunk=repmat(tmp_file_chunk,[max(unique(tmp_file_label)) 1]);
        cfg2.design = make_design_similarity(cfg2);
        cfg2.basic_checks.DoubleFilenameEntriesOk = 1;
    end
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
                
                if Trial_type==2
                    g2=0;
                    for Num_blocks=[(class-1)*(length(cfg2.files.name)./categories)+1:class*(length(cfg2.files.name)./categories)]
                        g2=g2+1;
                        tmp_beta_cat2=tmp_beta_cat2+niftiread(cfg2.files.name{Num_blocks,1}(1:end));
                    end
                    tmp_beta_cat2=tmp_beta_cat2./g2;
                    new_mean_file=fullfile([cfg2.files.name{Num_blocks,1}(1:end-13),'mean_c',sprintf('%0.3d', class),'_beta.nii']);
                    niftiwrite(tmp_beta_cat2,new_mean_file);
                    cfg2_tmp.files.name{class,1}=new_mean_file;
                    cfg2_tmp.files.chunk(class,1)=class;
                    cfg2_tmp.files.label(class,1)=class;
                    cfg2_tmp.files.descr{1,class}=cfg.files.descr{1,class};
                end
            end
            cfg_tmp.files.set=[];
            cfg_tmp.files.xclass=[];
            cfg.files=cfg_tmp.files;
            cfg.design.label=[1:categories]';
            cfg.design.train(categories+1:end)=[];
            cfg.design.test(categories+1:end)=[];
            
            if Trial_type==2
                cfg2_tmp.files.set=[];
                cfg2_tmp.files.xclass=[];
                cfg2.files=cfg2_tmp.files;
                cfg2.design.label=[1:categories]';
                cfg2.design.train(categories+1:end)=[];
                cfg2.design.test(categories+1:end)=[];
                
                cfg.files.name=vertcat(cfg.files.name,cfg2.files.name);
                cfg.files.chunk=vertcat(cfg.files.chunk,cfg2.files.chunk);
                cfg.files.label=vertcat(cfg.files.label,cfg2.files.label);
                cfg.design.label=vertcat(cfg.design.label,cfg.design.label);
                cfg.design.train=vertcat(cfg.design.train,cfg.design.train);
                cfg.design.test=vertcat(cfg.design.test,cfg.design.test);
                cfg.files.descr=horzcat(cfg.files.descr,cfg2.files.descr);
            end
        end
        
        if searchlighting==1
            cfg.analysis = 'searchlight';
            cfg.searchlight.unit = 'voxels';
            cfg.searchlight.radious = 100;
            cfg.searchlight.spherical = 1;
            cfg.plot_selected_voxels = 10;
            results = decoding(cfg);
            RSA_results=results;
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
                    cfg.files.mask = {[fullfile(Main_analysis_directory, 'Results_temp_playing'),'\ROIs\',csub,'\',Regions{region},'_roi.img']};                    
                    RSA_results{region}=decoding(cfg);
                    if Trial_type==2
                        RSA_results{region}.other.output{1, 1}=RSA_results{region}.other.output{1, 1}(categories+1:2*categories,1:categories);
                    end
                    close all;
                end
            end
        end
        
    elseif Decod_instead_of_Corr==1
        
        cfg_tmp=cfg;
        
        if Trial_type==2
            cfg_tmp2=cfg2;
        end
        cfg.results.output = {'accuracy','accuracy_minus_chance','balanced_accuracy','sensitivity','specificity','AUC','dprime','loglikelihood','corr','zcorr'}; %, 'SVM_pattern'};
        cfg.decoding.software = 'libsvm';
        cfg.decoding.method = 'classification';
        cfg.analysis = 'roi';
        cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification % parameters for libsvm (linear SV classification, cost = 1, no screen output)
        
        uniqe_labels=unique(cfg.files.label);
        combinations=nchoosek([1:length(uniqe_labels)],2);
        for comb=1:length(combinations)
            ind1=(cfg_tmp.files.label==combinations(comb,1));
            if Trial_type==1 || Trial_type==3 || Trial_type==4 || Trial_type==5 || Trial_type==6 || Trial_type==7 || Trial_type==8 || Trial_type==9 || Trial_type==10 || Trial_type==11 || Trial_type==12
                ind2=(cfg_tmp.files.label==combinations(comb,2));
            elseif Trial_type==2
                ind2=(cfg_tmp2.files.label==combinations(comb,2));
            end
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
            if Trial_type==1 || Trial_type==3 || Trial_type==4 || Trial_type==5 || Trial_type==6 || Trial_type==7 || Trial_type==8 || Trial_type==9 || Trial_type==10 || Trial_type==11 || Trial_type==12
                for n=1:length(ind2)
                    if ind2(n)==1
                        t=t+1;
                        name_tmp{t,1}=cfg_tmp.files.name{n};
                        chunk_tmp(t,1)=cfg_tmp.files.chunk(n);
                        label_tmp(t,1)=cfg_tmp.files.label(n);
                        desc_tmp{1,t}=cfg_tmp.files.descr{1,n};
                    end
                end
            elseif Trial_type==2
                for n=1:length(ind2)
                    if ind2(n)==1
                        t=t+1;
                        name_tmp{t,1}=cfg_tmp2.files.name{n};
                        chunk_tmp(t,1)=cfg_tmp2.files.chunk(n);
                        label_tmp(t,1)=cfg_tmp2.files.label(n);
                        desc_tmp{1,t}=cfg_tmp2.files.descr{1,n};
                    end
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
                        cfg_tmp_final.files.mask = {[fullfile(Main_analysis_directory, 'Results_temp_playing'),'\ROIs\',csub,'\',Regions{region},'_roi.img']};
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
    if searchlighting==1
        SL='SearchLight';
    else
        SL='ROIs';
    end
    if num_iterations>1
        Iter='Iter_';
    else
        Iter='';
    end
    if Trial_type==1
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==2
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name_scnd,'/',['RSA_',Iter,'Cor_Incor_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==3
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Cue_Color_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==4
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Stm_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==5
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Stm_clpsdCueCol_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==5
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Stm_Side_clpsdCueCol_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==6
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Stm_Side_clpsdCueCol_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==7
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Resp_',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==8
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Resp_acros_Cues',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==9
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Resp_opst_Colrs',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==10
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Stim_side_acrs_rules',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==11
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Resp_InOut_opst_Colrs',dec,SL,'.mat']),'RSA_results');
    elseif Trial_type==12
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_',Iter,'Resp_InOut_opst_Colrs2',dec,SL,'.mat']),'RSA_results');
    end
end
done=1;
end