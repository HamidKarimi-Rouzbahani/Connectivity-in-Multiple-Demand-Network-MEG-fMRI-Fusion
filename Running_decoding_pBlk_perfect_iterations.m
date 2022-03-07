clc;
clear all;
% close all;

tic
for subject=[1:26 28:31]
    
    Action=4;       % =1 for decoding of stim, rule and resp each two class
    % 2:16*16,   3:8*8    4:4*4 RDMs
    
    % 10 1 3 11: side, periphery, rule, resp
    
    Trial_type=10; % 1=correct (for RSA=colors*rules*stimuli; for decoding= regresors are periph, rules, resp ); 2=same errors (for RSA= the same as previous);
    % 3=swapped errors (only for decoding) % 4=correct (for decoding only= regresors are periph, side, rules, resp)
    
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
    
    
    Decod_instead_of_Corr=1; % Decoding rates or correlations in RSA? (only for RSA)
    Searchlighting=0;        % Searchlight or RoI? (for both rsa and decoding)
    
    
    Behavioural_directory='/group/woolgar-lab/projects/Hamid/Projects/MD/fMRI_Behavioural_data/';
    fmri_data_address='/group/woolgar-lab/projects/Hamid/Projects/MD/fMRI_Scans';
    Main_analysis_directory='/group/woolgar-lab/projects/Hamid/Projects/MD/Analyses/';
    
    
    % addpath(fullfile([Main_analysis_directory,'\spm12\']));
    
    %         Events_prepared=Step1_Event_preparation_f_pBlk(subject,Behavioural_directory,Main_analysis_directory);
    %
    % Files_converted=Step2_DICOM_to_NIFTI_f(subject,fmri_data_address,Main_analysis_directory);
    %
    % Preprocessing_done=Step3_Preprocessing_running_f(subject,Behavioural_directory,Main_analysis_directory);
    
    %     Regression_definition_done=Step4_1st_level_regression_new_f_pBlk(subject,Main_analysis_directory,Behavioural_directory,Action,Trial_type);
    %     
    %     Rgression_B_estimate_done=Step5_Beta_Weights_Estimation_f_pBlk(subject,Main_analysis_directory,Action,Trial_type);
    %
    %     %     MD_regions_warping_done=Step6_Unwarp_MD_regions_f_pBlk(subject,Main_analysis_directory,Action,Trial_type)
    %
    %     %     My_specific_regions_warping_done=Step6_5_Unwarp_my_regions_f_pBlk(subject,Main_analysis_directory,Action,Trial_type)
    %
    %     %     Visual_regions_warping_done=Step7_Unwarp_visual_regions_f_pBlk(subject,Main_analysis_directory,Action,Trial_type)
    %
    Decoding_done=Step8_Decoding_f_pBlk_iterations(subject,Main_analysis_directory,Action,Trial_type,Decod_instead_of_Corr,Searchlighting);
    %
    %                 Decoding_done=Step9_Decoding_specific_ROIs_f_pBlk(subject,Main_analysis_directory,Action,Trial_type,Decod_instead_of_Corr,Searchlighting);
    [subject]
end
toc