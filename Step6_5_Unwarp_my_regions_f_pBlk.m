function [done]=Step6_5_Unwarp_my_regions_f_pBlk(subject,Main_analysis_directory,Action,Trial_type)

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
end

% AWoolgar
% do unwarp BA4 region and save to subject ROI folder

% update Sept 2009: will skip subject if unwarped MD of that name already
% exists

% clear

%rootdir = '/imaging/dm06/TMS-fMRI';

rootdir = fullfile([Main_analysis_directory,'\Results_temp_playing']);
dataroot = fullfile(rootdir); %where the subject directories are
maskdir = fullfile(rootdir);



% add marsbar to the path
addpath(fullfile([Main_analysis_directory,'\spm12\toolbox\marsbar']));
addpath(fullfile([Main_analysis_directory,'\spm12']));

%spm fmri

% load data
% resmat =  fullfile(rootdir,'behav', 'first7.mat'); %the resmat with everybody's data in it
% load (resmat); % struc = res
% subs = res.subs;
% nsubs = length(subs);


%roidir ='/imaging/dm06/TMS-fMRI/tools/brodmann_areas/';
roidir=fullfile([Main_analysis_directory,'\My_ROIs']);

rois = {'Above_chance_decoding_rule_clusters_1_31_39_roi',...
    'Above_chance_decoding_rule_clusters_29_18_33_roi',...
    'Above_chance_decoding_stim_clusters_7_-77_17_roi',...
    'Above_chance_decoding_stim_clusters_46_17_26_roi'};


skipem = [1];

for s = subject%:nsubs %loop over Ss
    csub = sprintf('%s%0.3d','S', s);
    
    disp(csub); %['rd' sprintf('%0.3d', s)]; %alls(s,:);
    
    outdir = fullfile(rootdir, 'ROIs' , csub);
    
    for r=1:length(rois)
        
        % name ROI to be created
        outfile = char(fullfile(outdir, strcat(rois(r), '.img')));
        %  outfile = fullfile(outdir, [roitype 'MD' num2str(r) '.img'])
        
        % outfile = fullfile(outdir, ['BA1718_' num2str(r) '.img']);
        
        % only do this loop if the unwarped ROI does not already exist!
        % if exist(outfile); disp('skipping'); continue; end
        
        %make output dir if it doesn't exist
        if exist(fullfile(outdir)) ~=7; mkdir(fullfile(outdir)); end
        
        % -- prepare outimage --
        baseimg = fullfile(maskdir, csub,Beta_folder_name,'Block001', 'mask.nii');
        
        MMmask = spm_vol(baseimg);
        img = MMmask;
        img.fname = outfile;
        img.descrip = 'Unwarped ROI for this subject'; %
        vals = NaN(img.dim(1),img.dim(2),img.dim(3));
        
        % -- unnormalise voxels --
        snmatfile = dir(fullfile(dataroot, csub, 'structurals', '*seg_sn.mat')); %normalisation params file in structurals dir
        %check only found one option!
        if size(snmatfile,1) == 1 %if it just found 1
            snmat = fullfile(dataroot, csub, 'structurals', snmatfile.name); %normalisation params
        else
            error('size(snmatfile,1)~=1, maybe you haven''t calculated normalisation params for this subject, or maybe there are two or more files ending in _sn.mat in this subject''s structurals directory? OR the folder name is wrongly defined')
        end
        croi = maroi(fullfile(roidir,(rois{r}))); %make marsy object from current roi
        vox = realpts(croi,native_space(croi)); % points in mm
        
        % do deformation of ROI for current subject:
        vox = round(spm_get_orig_coord(vox',snmat, baseimg)); %in voxels
        [uvox vi] = unique(cellstr(num2str(vox')));
        vox = vox(:,sort(vi));
        vox = vox';
        
        % - cross reference with whole brain mask -
        mY = spm_sample_vol(MMmask,vox(1,:),vox(2,:),vox(3,:),1);
        vox = vox(:,find(mY==1));
        
        % - place 1s in those co-ords in vals -
        for v = 1:size(vox,2)
            vals(vox(1,v),vox(2,v),vox(3,v)) = 1;
        end
        
        % - write image -
        spm_write_vol(img,vals);
    end
end

done=1;
end