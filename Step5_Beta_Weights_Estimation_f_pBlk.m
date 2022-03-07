function [done]=Step5_Beta_Weights_Estimation_f_pBlk(subject,Main_analysis_directory,Action,Trial_type)

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

csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
csub = csub(1:4);

folders=dir(fullfile(Main_analysis_directory, 'Results_temp_playing', csub,Beta_folder_name));
c=0;
for blk_counter=1:size(folders,1)
    if strncmp(folders(blk_counter).name,'Bl',2)
        c=c+1;
    end
end
for blk_counter=1:c
    
    Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {[fullfile(Main_analysis_directory, 'Results_temp_playing', csub,'/',Beta_folder_name,num2str(blk_counter,'/Block%03d')),'\SPM.mat']};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    save(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,num2str(blk_counter,'/Block%03d'),'/',['Level1_Estimate_per_block_Subj_',num2str(subject),'.mat']),'matlabbatch');
    
    load([fullfile(Main_analysis_directory, 'Results_temp_playing', csub,'/',Beta_folder_name,num2str(blk_counter,'/Block%03d')),'\SPM.mat']);
    tmp=0;
    for i=1:length(SPM.Sess.U)
        tmp=tmp+length(SPM.Sess.U(i).ons);
    end
    if tmp>0
        spm_jobman('run',matlabbatch);
    end
end

done=1;
end