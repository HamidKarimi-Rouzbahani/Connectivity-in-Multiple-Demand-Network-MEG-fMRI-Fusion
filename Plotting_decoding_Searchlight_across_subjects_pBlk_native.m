clc;
% close all;
clear all;
Main_analysis_directory='D:/Hamid/Postdoc/MD/Analyses';

Action=1;       % =1 for decoding of stim, rule and resp each two class
%2:16*16,   3:8*8    4:4*4 RDMs
Trial_type=2; % 1=correct; 2=all errors


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
        
infos={'S','R'};
f=0;
for info=[1:2]
    for subject=subjects
        f=f+1;
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);        
        load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_',SL,'_',infos{info},'.mat']));
        beta_file_temp=niftiread(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/Block001/beta_0001.nii'));
        image_info = niftiinfo(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/Block001/beta_0001.nii'));
ccc
        beta_file_temp(Decoding_results.mask_index)=Decoding_results.accuracy_minus_chance.output;
        if Trial_type==2
            beta_file_temp=-beta_file_temp;
        end
        new_file=fullfile([Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/Decoding_searchlight_image_',infos{info},'.nii']);
%         imagesc(beta_file_temp(:,:,15)); colorbar
%         pause (2)
        niftiwrite(beta_file_temp,new_file,image_info);
        norm_params = dir(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/','structurals','*seg_sn.mat')); %normalisation params file in structurals dir
        norm_params = fullfile(Directory_for_working,num2str(csub,'S%03d'),'/','structurals',norm_params.name); 
        spm_write_sn(new_file,norm_params,options);
        [info subject]
    end
end
ccc
%% run the group-level subjects decoding searchlight
clc;
infos={'S','R'};
info=2; % 1=Stim; 2=Rule
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.dir = {['D:\Hamid\Postdoc\MD\Analyses\Results_temp_playing\Group_level_analyses\',Beta_folder_name,infos{info}]};

s=0;
for subject=subjects
    s=s+1;
    csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
    csub = csub(1:4);
    matlabbatch{1, 1}.spm.stats.factorial_design.des.t1.scans{s,1}=fullfile([Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/wDecoding_searchlight_image_',infos{info},'.nii']);
end
mkdir(fullfile([Directory_for_working,'/Group_level_analyses/',Beta_folder_name,infos{info},'/']))
save(fullfile([Directory_for_working,'/Group_level_analyses/',Beta_folder_name,infos{info},'/',infos{info},'.mat']),'matlabbatch');
spm_jobman('run',matlabbatch);

clearvars matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile([Directory_for_working,'/Group_level_analyses/',Beta_folder_name,infos{info},'/','SPM.mat'])};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
spm fmri
