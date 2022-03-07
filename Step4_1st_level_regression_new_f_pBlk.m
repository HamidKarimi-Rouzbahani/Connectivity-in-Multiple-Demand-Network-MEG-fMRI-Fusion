function [done]=Step4_1st_level_regression_new_f_pBlk(subject,Main_analysis_directory,Behavioural_directory,Action,Trial_type)



if Action==1
    Beta_folder_name='Correct_trials_decoding_design_pBlk_native';
elseif Action==2
    Beta_folder_name='Correct_trials_RSA_design_16_pBlk_native';
elseif Action==3
    Beta_folder_name='Correct_trials_RSA_design_8_pBlk_native';
elseif Action==4
    Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
end
if Trial_type==2
        Beta_folder_name=[Beta_folder_name,'_AllErr'];
elseif Trial_type==3
    if Action==1
        Beta_folder_name=[Beta_folder_name,'_OpsErr'];
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

native_space=1;
Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];
addpath(fullfile([Main_analysis_directory,'\spm12\']));

runss=dir([Behavioural_directory,'/test_res_p',num2str(subject),'_r*.mat']);

for i=[1:size(runss,1)]
    epi{i,1} = ['run',num2str(i)];
end
struc = 'structurals';

csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
csub = csub(1:4);
mkdir(fullfile([Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/']));

if Action==1
    if Trial_type==1    
        desired_condition_names={'S1','S2','R1','R2','P1','P2'};
    elseif Trial_type==4
        desired_condition_names={'S1','S2','R1','R2','P1','P2','L1','L2'}; % L refers to stimulus sides (location)
    end
    no_of_first_blk_conds=2;
    % Coding of names
    % First char S:stimulus; R:Rule; C:Rule color; P:Response
    % Second char: conition
    % Third char C:correct; S:stimulu error; R:rule error; B:both error
    
elseif Action==2
    desired_condition_names={'C1S1','C1S2','C1S3','C1S4',...
        'C2S1','C2S2','C2S3','C2S4',...
        'C3S1','C3S2','C3S3','C3S4',...
        'C4S1','C4S2','C4S3','C4S4'};
    no_of_first_blk_conds=16;
    
elseif Action==3
    desired_condition_names={'R1S1','R1S2','R1S3','R1S4',...
        'R2S1','R2S2','R2S3','R2S4'};
    no_of_first_blk_conds=8;
    
elseif Action==4
    if Trial_type==1 
        desired_condition_names={'R1S1','R1S2','R2S1','R2S2'};
    elseif Trial_type==3
        desired_condition_names={'R1C1','R1C2','R2C1','R2C2'};
    elseif Trial_type==4 || Trial_type==5 || Trial_type==6
        desired_condition_names={'S1','S2','S3','S4'};
    elseif Trial_type==7
        desired_condition_names={'R1','R2','R3','R4'};
    elseif Trial_type==8
        desired_condition_names={'R1P1','R1P2','R2P1','R2P2'};
    elseif Trial_type==9
        desired_condition_names={'C1P1','C1P2','C2P1','C2P2'};        
    elseif Trial_type==10
        desired_condition_names={'C1S1','C1S2','C2S1','C2S2'};
    elseif Trial_type==11
        desired_condition_names={'C1P1','C1P2','C2P1','C2P2'};
    elseif Trial_type==12
        desired_condition_names={'C1P1','C1P2','C2P1','C2P2'};        
    end
    no_of_first_blk_conds=4;
end

for i=1:size(desired_condition_names,2)
    if Trial_type==1
        desired_condition_names{1,i}(end+1)='C';
    elseif Trial_type==2
        if Action==1
            desired_condition_names{1,i}(end+1)=desired_condition_names{1,i}(1);
        elseif Action>1
            desired_condition_names{1,i}(end+1)='E';
        end
    elseif Trial_type==3 
        if Action==1
            if strcmp(desired_condition_names{1,i}(1),'S')
                desired_condition_names{1,i}(end+1)='R';
            elseif strcmp(desired_condition_names{1,i}(1),'R')
                desired_condition_names{1,i}(end+1)='S';
            elseif strcmp(desired_condition_names{1,i}(1),'P')
                desired_condition_names{1,i}(end+1)='E';
            end
        elseif Action>1
            desired_condition_names{1,i}(end+1)='C';
        end
    elseif Trial_type==4 || Trial_type==5 || Trial_type==6 || Trial_type==7 || Trial_type==8 || Trial_type==9 || Trial_type==10 || Trial_type==11 || Trial_type==12
        desired_condition_names{1,i}(end+1)='C';
    end
end

blk_counter=0;
for run=1:length(epi)
    clear my_epi %no hangovers from prev subs
    if exist(fullfile(Main_analysis_directory, 'Results_temp_playing', csub, epi{run})) %only do this for folders that exist - not all people have all runs
        my_epi{run} = fullfile(Main_analysis_directory, 'Results_temp_playing', csub, epi{run});
    end
    my_struc = fullfile(Main_analysis_directory, 'Results_temp_playing', csub, struc);
    
    if native_space==1
        files = get_files(char(my_epi{run}),'arf*.nii');
    else
        files = get_files(char(my_epi{run}),'s4warf*.nii');
    end
    
    load(fullfile(Directory_for_working,num2str(csub,'S%03d'),['Behav_','Subj_',num2str(subject),'_run_',num2str(run),'_per_blk.mat']));
    % rows 1:4==> stim err, both err, rule err, correct
    % colums 1:4==> conditions
    
    
    stim_times=[];
    stim_dures=[];
    for i=1:4
        for j=1:2
            if ~isempty(Stim{i,j})
                stim_times=vertcat(stim_times,Stim{i,j});
                stim_dures=vertcat(stim_dures,Stim_rt{i,j});
            end
        end
    end
    [mx,no]=max(stim_times);
    last_slice_first_block=ceil((mx+stim_dures(no)+21)./2);
    
    for blk=1:2
        matlabbatch{1, 1}.spm.stats.fmri_spec.timing.units='secs';
        matlabbatch{1, 1}.spm.stats.fmri_spec.timing.RT=2;
        matlabbatch{1, 1}.spm.stats.fmri_spec.timing.fmri_t=34;
        matlabbatch{1, 1}.spm.stats.fmri_spec.timing.fmri_t0=1;
        blk_counter=blk_counter+1;
        if blk==1
            image_nums=1:last_slice_first_block;
            time_deduction=0;
        else
            image_nums=last_slice_first_block+1:size(files,1);
            time_deduction=last_slice_first_block*2;
        end
        in=0;
        for image_num=image_nums
            in=in+1;
            matlabbatch{1,1}.spm.stats.fmri_spec.sess(1).scans{in,1}=files(image_num,:);
        end
        
        if Action==1
            for condit=[1:length(desired_condition_names)]
                
                if strcmp(desired_condition_names{condit}(end),'C')
                    rows=4;
                elseif strcmp(desired_condition_names{condit}(end),'E')
                    rows=1:3;
                elseif strcmp(desired_condition_names{condit}(end),'S')
                    rows=1;
                elseif strcmp(desired_condition_names{condit}(end),'R')
                    rows=3;
                elseif strcmp(desired_condition_names{condit}(end),'B')
                    rows=2;
                end
                
                onsets=[];
                durations=[];
                for row=rows
                    if strcmp(desired_condition_names{condit}(1),'S')
                        onsets=vertcat(onsets,Stim{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*no_of_first_blk_conds]}-time_deduction);
                        durations=vertcat(durations,Stim_rt{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*no_of_first_blk_conds]});
                    elseif strcmp(desired_condition_names{condit}(1),'R')
                        onsets=vertcat(onsets,Rule{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*2]}-time_deduction);
                        durations=vertcat(durations,Rule_rt{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*2]});
                    elseif strcmp(desired_condition_names{condit}(1),'C')
                        onsets=vertcat(onsets,Rule_color{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*4]}-time_deduction);
                        durations=vertcat(durations,Rule_color_rt{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*4]});
                    elseif strcmp(desired_condition_names{condit}(1),'P')
                        onsets=vertcat(onsets,Respdec{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*no_of_first_blk_conds]}-time_deduction);
                        durations=vertcat(durations,Resp_rtdec{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*no_of_first_blk_conds]});
                    elseif strcmp(desired_condition_names{condit}(1),'L')
                        onsets=vertcat(onsets,StimTside{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*no_of_first_blk_conds]}-time_deduction);
                        durations=vertcat(durations,Stim_rtTside{row,[str2num(desired_condition_names{condit}(2))+(blk-1)*no_of_first_blk_conds]});
                    end
                end
                
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).name=desired_condition_names{condit};
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).onset=onsets;
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).duration=durations;
                
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).tmod=0;
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).pmod=struct([]);
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).orth=1;
                if isempty(matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).onset)
                    missing_regressors(condit)=1;
                else
                    missing_regressors(condit)=0;
                end
            end
            
        elseif Action>1
            for condit=[1:length(desired_condition_names)]
                
                if strcmp(desired_condition_names{condit}(end),'C')
                    rows=4;
                elseif strcmp(desired_condition_names{condit}(end),'E')
                    rows=1:3;
                elseif strcmp(desired_condition_names{condit}(end),'S')
                    rows=1;
                elseif strcmp(desired_condition_names{condit}(end),'R')
                    rows=3;
                elseif strcmp(desired_condition_names{condit}(end),'B')
                    rows=2;
                end
                
                if Trial_type==1
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Stim_Col_4;
                        Data_matrix_rt=Stim_Col_rt_4;
                    end
                elseif Trial_type==3
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Rule_color;
                        Data_matrix_rt=Rule_color_rt;
                    end
                elseif Trial_type==4
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Stim_4_stims_2_Blk;
                        Data_matrix_rt=Stim_rt_4_stims_2_Blk;
                    end
                elseif Trial_type==5
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Stim_Col_Cues_Colpsd_4;
                        Data_matrix_rt=Stim_Col_Cues_Colpsd_rt_4;
                    end
                elseif Trial_type==6
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Stim_side_Col_Cues_Colpsd_4;
                        Data_matrix_rt=Stim_side_Col_Cues_Colpsd_rt_4;
                    end
                elseif Trial_type==7
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Resp;
                        Data_matrix_rt=Resp_rt;
                    end
                elseif Trial_type==8
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Resp_Col_4;
                        Data_matrix_rt=Resp_Col_rt_4;
                    end
                elseif Trial_type==9
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Resp_Col_Cues_Colpsd_4;
                        Data_matrix_rt=Resp_Col_Cues_Colpsd_rt_4;
                    end
                elseif Trial_type==10
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Stim_side_Col_4;
                        Data_matrix_rt=Stim_side_Col_rt_4;
                    end                   
                elseif Trial_type==11
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Resp_InOut_Cues_Colpsd_4;
                        Data_matrix_rt=Resp_InOut_Cues_Colpsd_rt_4;
                    end
                elseif Trial_type==12
                    if Action==2
                        Data_matrix=Stim_Col_16;
                        Data_matrix_rt=Stim_Col_rt_16;
                    elseif Action==3
                        Data_matrix=Stim_Col_8;
                        Data_matrix_rt=Stim_Col_rt_8;
                    elseif Action==4
                        Data_matrix=Resp_InOut_Cues_Colpsd2_4;
                        Data_matrix_rt=Resp_InOut_Cues_Colpsd2_rt_4;
                    end                    
                end
                onsets=[];
                durations=[];
                for row=rows
                    onsets=vertcat(onsets,Data_matrix{row,[condit+(blk-1)*no_of_first_blk_conds]}-time_deduction);
                    durations=vertcat(durations,Data_matrix_rt{row,[condit+(blk-1)*no_of_first_blk_conds]});
                end
                
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).name=desired_condition_names{condit};
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).onset=onsets;
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).duration=durations;
                
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).tmod=0;
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).pmod=struct([]);
                matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).orth=1;
                if isempty(matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).cond(condit).onset)
                    missing_regressors(condit)=1;
                else
                    missing_regressors(condit)=0;
                end
            end
        end
        
        matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).multi={''};
        matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).regress=struct([]);
        matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).multi_reg={''};
        matlabbatch{1, 1}.spm.stats.fmri_spec.sess(1).hpf=128;
        matlabbatch{1, 1}.spm.stats.fmri_spec.fact=struct([]);
        matlabbatch{1, 1}.spm.stats.fmri_spec.bases.hrf.derivs=[0 0];
        matlabbatch{1, 1}.spm.stats.fmri_spec.volt=1;
        matlabbatch{1, 1}.spm.stats.fmri_spec.global='None';
        matlabbatch{1, 1}.spm.stats.fmri_spec.mthresh=0.8;
        matlabbatch{1, 1}.spm.stats.fmri_spec.mask=cell(1,1);
        matlabbatch{1, 1}.spm.stats.fmri_spec.cvi='AR(1)';
        matlabbatch{1, 1}.spm.stats.fmri_spec.dir{1,1}=fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name,num2str(blk_counter,'/Block%03d'));
        mkdir(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,num2str(blk_counter,'/Block%03d')))
        save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,num2str(blk_counter,'/Block%03d'),'/',['Level1_Regress_per_blk_Subj_',num2str(subject),'.mat']),'matlabbatch','missing_regressors');
        spm_jobman('run',matlabbatch);
        clearvars matlabbatch
        [blk_counter]
    end
end

done=1;


end