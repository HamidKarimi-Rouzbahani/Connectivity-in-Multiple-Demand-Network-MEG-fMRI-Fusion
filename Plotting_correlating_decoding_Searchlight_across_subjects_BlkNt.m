clc;
% close all;
clear all;
Main_analysis_directory='D:/Hamid/Postdoc/MD/Analyses';

Action=1;       % =1 for decoding of stim, rule and resp each two class
%2:16*16,   3:8*8    4:4*4 RDMs
Trial_type=1; % 1=correct; 2=all errors


Decod_instead_of_Corr=1;
searchlighting=1;

subjects=[1:26 28:31];

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


if Action==1
    Beta_folder_name='Correct_trials_decoding_design_pBlk_native';
elseif Action==2
    Beta_folder_name='Correct_trials_RSA_design_16_pBlk_native';
    MEG_RDM=RDM_correct_16_all;
    NoC=16;
elseif Action==3
    Beta_folder_name='Correct_trials_RSA_design_8_pBlk_native';
    MEG_RDM=RDM_correct_8_all;
    NoC=8;
elseif Action==4
    Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
    MEG_RDM=RDM_correct_4_all;
    NoC=4;
end
if Trial_type==2
    Beta_folder_name=[Beta_folder_name,'_AllErr'];
elseif Trial_type==3
    Beta_folder_name=[Beta_folder_name,'_OpsErr'];
end


Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];

options.bb = [  -78  -112   -50; 78    76    85];   % = bounding box = the range of
options.vox = [2 2 2];                              % size of voxels to use in the normalised images
options.interp = 1;                                 % interpolation method
options.wrap = [0 0 0];                             % wrap edges?
options.preserve = 0;                               % preserve voxel concentrations? (mainly for VBM -

load('behavioural_results.mat')
infos={'S','R'};
f=0;
Warped_decoding_file=nan(79,95,68,31);
beh=2; % 1= accuracy; 2=RT
if beh==1
    Thing=Stims_no/0.4;
else
    Thing=Stims_rt;
end
trial_outcome=4;
for info=[1:2]
    for subject=subjects
        f=f+1;
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        Warped_decoding_file(:,:,:,subject)=niftiread(fullfile([Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/wDecoding_searchlight_image_',infos{info},'.nii']));
        image_info{subject}=niftiinfo(fullfile([Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/wDecoding_searchlight_image_',infos{info},'.nii']));
        %         imagesc(Warped_decoding_file(:,:,15,subject)); colorbar
        %         pause (2)
        [info subject]
    end
    for x=1:size(Warped_decoding_file,1)
        for y=1:size(Warped_decoding_file,2)
            for z=1:size(Warped_decoding_file,3)
                beta_file_temp(x,y,z)=corr(squeeze(Warped_decoding_file(x,y,z,subjects)),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,subjects)),3),2)));
            end
        end
    end
%     imagesc(beta_file_temp(:,:,15)); colorbar
    
    mkdir(fullfile([Directory_for_working,'/Group_level_analyses/',Beta_folder_name,'_across/']))
    if beh==1
        new_file=fullfile([Directory_for_working,'/Group_level_analyses/',Beta_folder_name,'_across/wBeh_Acc_corr_searchlight_image_',infos{info},'.nii']);
    elseif beh==2
        new_file=fullfile([Directory_for_working,'/Group_level_analyses/',Beta_folder_name,'_across/wBeh_RT_corr_searchlight_image_',infos{info},'.nii']);
    end
    niftiwrite(single(beta_file_temp),new_file,image_info{subject});
end
spm fmri
