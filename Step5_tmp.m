

subject=1;

Action=4;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs

Trial_type=12; % 1=correct (for RSA=colors*rules*stimuli); 2=same errors (for RSA= the same as previous);
% 3=swapped errors (only for decoding) % 4=correct (for decoding=colors*rules*stimuli*stim side)

% 3=rule coding across cue colors (only for RSA: 1 to 4 cue colors collapsed across all stimuli)
% 4=stimulus coding  (only for RSA: 1 to 4 stimuli collapsed across rules and cues):
% 5=peripheral stimulus coding (only for RSA: collapsed across stimulus side and cue colors from different rules)
% 6=side stimulus coding (only for RSA: collapsed across cue colors from different rules)
% 7=Response coding (responses 1-4)
% 8=Response coding (with rules included: columns are rule1resp1 rule1resp2 rule2resp1 rule2resp2)
% 9=Response coding (with rules included: columns are color1resp1 color1resp2 color2resp1 color2resp2)
% 10=Stim side coding (with rules included: columns are color1stim1 color1stim2 color2stim1 color2stim2)
% 11=Response coding Inner vs Outer and not left and right(with rules included: columns are color1resp1 color1resp2 color2resp1 color2resp2)
% 12=Response coding Inner vs Outer and not left and right(same as 11 but collapsed different colors with rules included: columns are color1resp1 color1resp2 color2resp1 color2resp2)

Decod_instead_of_Corr=1;
searchlighting=0;


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