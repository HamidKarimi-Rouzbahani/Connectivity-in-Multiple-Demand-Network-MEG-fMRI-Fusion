function [done]=Step3_Preprocessing_running_f(subject,Behavioural_directory,Main_analysis_directory)
% HKR March 2020, do preprocessing

% from AW dec 2018, do preprocessing

% Originally adapted from R Thompson's example script

runss=dir([Behavioural_directory,'/test_res_p',num2str(subject),'_r*.mat']);

for i=[1:size(runss,1)]
    epi{i} = ['run',num2str(i)];
end
rootdir = Main_analysis_directory;

struc = 'structurals';

% latest version of SPM:
addpath(fullfile([Main_analysis_directory,'\spm12']));
addpath(fullfile([Main_analysis_directory,'\spm12\toolbox\spmtools-tsdiffana-code']));
addpath(fullfile([Main_analysis_directory,'\spm12\toolbox\tsdiffana']));

% addpath('F:\Postdoc\MD_project\fMRI\Data\Analyses'); %for compatibility with old spm_imcalc_ui

%spm fmri;

% n_dummy_scans = 0;

% % subject loop
% resdir =  fullfile(rootdir,'scanner_log_files');
% resmat = fullfile(resdir,'first43.mat');
% load(resmat);
% subs = res.subs;
% failedsubs = '';

overwrite = 1; %overwrite previous rather than skipping steps

csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
disp(csub);
csub = csub(1:4);
disp(csub);
addpath([Main_analysis_directory,'\Results_temp_playing\',csub])

%try

% Step 1 - set up dirs
clear my_epi %no hangovers from prev subs
for i = 1:length(epi)
    if exist(fullfile(rootdir, 'Results_temp_playing', csub, epi{i})); %only do this for folders that exist - not all people have all runs
        my_epi{i} = fullfile(rootdir, 'Results_temp_playing', csub, epi{i});
    end
end
my_struc = fullfile(rootdir, 'Results_temp_playing', csub, struc);

% name of print file that will be used to save various graphs etc generated
% as part of the preprocssing
% ------------------------------------------------------------------------
if ~exist(fullfile(rootdir, 'preprocessing')); mkdir(fullfile(rootdir, 'preprocessing')); end
printfile = fullfile(rootdir, 'preprocessing',['spm_preprocessing_' csub '.ps']);
pdffile = fullfile(rootdir, 'preprocessing',['spm_preprocessing' csub '.pdf']);

% =========================================================================
% =========================================================================
% Step  2 - Use tsdiffana to have a look at the EPI data and check for
% obvious artefacts etc
% =========================================================================

% if output file doesn't already exist...
if ~exist(fullfile(my_epi{1},'raw_mean_img.nii')) || overwrite==1;
    
    %disp('running tsdiffana');
    
    %get files for all sessions at once
    files = get_files(char(my_epi),'f*.nii');
%     for i=1:length(my_epi)
%         tsdiffana(get_files(char(my_epi{i}),'f*.nii')); % pass the file names into tsdiffana % can't
%         pause(10);
%     end
    %make tsdiffana work at MACCS yet - recursion limit error.
    %spm_print(printfile);
    
    % use SPM imcalc (a very useful function that allows you to perform
    % calculations on sets of images) to generate a mean image
    spm_imcalc_ui(files,fullfile(my_epi{1},'raw_mean_img.nii'),'mean(X)',{1;0;16;1});
    %        spm_imcalc_ui(files,fullfile(my_epi{1},'raw_mean_img.nii'),'var(X)',{1;0;16;1});
else
    disp('skipping step 2 as mean image exists');
    files = get_files(char(my_epi),'f*.nii');
end

% =========================================================================

% =========================================================================
% Step  3 - Realign the EPI data so its all in the same space
% =========================================================================

%if first re-aligned file does not already exist
[path name ext] = fileparts(files(1,:));
if ~exist(fullfile(path, ['r' name ext])) || overwrite==1;
    
    disp('running realign');
    
    % set up some options for the realignment procedure:
    %---------------------------------------------------
    options.quality = 0.9000; % quality of results (0 = worse + quicker, 1 =
    % better + slower)
    options.weight = 0; % use a weighting image? Weighting can be specified if
    % you
    options.interp = 2; % method of interpolation
    options.wrap = [0 0 0];
    options.sep = 4; % precision with which images are sampled in order to
    % determine the realignment parameters
    options.fwhm = 5; % how much to smooth data for the purpose of realigning
    % (the images themselves won't be smoothed)
    options.rtm = 0; % if 0, images will be realigned to the first image.
    %If 1, images will be realigned to first image, then a mean image will be
    % calculated, then all images will be aligned to that mean image
    
    
    % Do the realignment --------------------------------
    %----------------------------------------------------
    % now pass the names of all the files to the realignment routine - this
    % will realign the data from all sessions into the same space. The way
    % we've done it, it will realign all the images to the very first image.
    
    spm_realign(files,options);
    
    
    % print the realignment parameter graphs to the output file:
    spm_print(printfile)
    %spm_print(pdffile)
    
    
    % Now reslice the images ----------------------------
    %----------------------------------------------------
    
    % At this point, the realignment is saved only in the headers of the
    % 'f*.nii' files. Now we're going to write out a set of new files that
    % incorporate this transformation.
    
    options.mask = 1; % mask out voxels that are outside the brain in any of
    % the images (i.e. only include voxels that are inside the brain in all images)
    options.mean = 1; % create a mean image
    options.interp = 4; % interpolation method
    options.which = 2; % which images to reslice - 2 = all
    
    
    spm_reslice(files,options);
else
    disp('skipping step 3 realign and reslice as first resliced image exists')
end

% =========================================================================
%% Step  4 - Slice Timing: Correct differences in slice acquisition times
% output = arf*.nii
% See https://github.com/automaticanalysis/automaticanalysis for code

% From Denise code - 'AttDec_step2_preprocessing.m'
%if first slice time corrected file does not already exist

[path name ext] = fileparts(files(1,:));
if ~exist(fullfile(path, ['ar' name ext])) || overwrite==1;
    
    %      STcorr_files = dir(fullfile(data_dir,'Results_temp_playing',num2str(sub_num,'S%03d'),runfolder_name{1,1},'ar*.nii')); % find all files
    %      if ~size(STcorr_files,1) > 1 || overwrite==1;
    
    % Read in the headers
    %load(fullfile(rootdir,'Results_temp_playing',['S', num2str(csub,'S%03d')],['S',num2str(csub,'S%03d'),'_Header_Files.mat']));
    load(fullfile(rootdir,'Results_temp_playing',num2str(csub,'S%03d'),[num2str(csub,'S%03d'),'_Header_Files.mat']));
    
    runfolder_name = dir(fullfile(rootdir,'Results_temp_playing',num2str(csub,'S%03d'),'run*'));
    runfolder_name = {runfolder_name.name};
    
    % CS 19  Get the repetition time (in S) - header file units are ms
    TR = h(1).hdrs{1}.RepetitionTime/1000;
    
    % There seems to be a bit of a spread in the slicetimes we're finding.
    % So we can't just pick the times from the first scan as a reference.
    % Instead, we'll get the slice times from all of the scans, and pick the mode.
    
    
    % Loop over runs
    for run_num = 1:length(epi)
        % Get the file names for this run
        %              files2 = [];
        %              tmp = dir(fullfile(rootdir,'Results_temp_playing',['S' num2str(csub,'S%03d')],runfolder_name{1,run_num},'rf*.nii')); % find all files
        %              tmp = [repmat([fullfile(rootdir,'Results_temp_playing',['S' num2str(csub,'S%03d')],runfolder_name{1,run_num}),filesep],size(tmp,1),1),char(tmp.name)]; % build the full path name for these files
        %              files2 = char(files2,tmp);
        %              % Remove empty file name
        %              files2 = files2(~all(files2'==' ')',:);
        
        rf_files = get_files(char(my_epi{run_num}),'rf*.nii');
        
        % Change the directory
        cd(fullfile(rootdir,'Results_temp_playing',num2str(csub,'S%03d'),runfolder_name{1,run_num}));
        
        % Loop over scans
        slicetimes_mat = [];
        for scan_num = 1:size(h(run_num).hdrs,2)
            
            str = h(run_num).hdrs{scan_num}.CSAImageHeaderInfo;
            idx = strcmp({str.name},'MosaicRefAcqTimes');
            
            % Get the slice times
            s_t = nan(str(idx).nitems,1);
            for j=1:str(idx).nitems
                s_t(j,1) = str2double(str(idx).item(j).val);
            end
            
            % Remove NaN values
            slicetimes = s_t(isfinite(s_t));
            
            % Add to the matrix
            slicetimes_mat = [slicetimes_mat,slicetimes];
        end
        
        % Now we can pick the mode of the slicetimes
        slicetimes = nan(size(slicetimes_mat,1),1);
        for slice_num = 1:size(slicetimes_mat,1)
            % Get the number of unique slicetime
            pos_times = unique(slicetimes_mat(slice_num,:));
            % See how often they occur
            n_occ = zeros(1,size(pos_times,2));
            for i = 1:size(pos_times,2)
                n_occ(1,i) = sum(slicetimes_mat(slice_num,:)==pos_times(1,i));
            end
            % Pick the slicetimes with the most occurences
            [~,idx] = max(n_occ);
            slicetimes(slice_num,1) = pos_times(1,idx);
        end
        
        
        
        %            % Now we can do the slice timing correction
        %spm_12_slice_timing(files,slicetimes,0,[0,TR])
        spm_slice_timing(rf_files,slicetimes,0,[0,TR]) %AW 14/11/19: 3rd argument is time in ms at which ref slice occurred, ie 0 means realigning everything to the start of the TR
        
    end
else
    disp('skipping step 4 slice timing as output file exists');
end


%     % =========================================================================
%     % Step  4 - Slice Timing
%     % =========================================================================
%
%
%     %if first slice time corrected file does not already exist
%     [path name ext] = fileparts(files(1,:));
%     if ~exist(fullfile(path, ['ar' name ext])) || overwrite==1;
%
%         disp('running slice timing');
%
%         %options:
%         sliceorder = [1:2:35 2:2:34]; % order in which slices are acquired - UReading "interleaved" - depends on n slices - see https://practicalfmri.blogspot.com.au/2012/07/siemens-slice-ordering.html or pdf in resources folder
%
%         refslice = 1; % which slice to align the timings to, in this case the first
%         % slice to be acquired
%         timings(1) = (2960-1000)/1000/35; % time to acquire one slice - in this case the TR (2
%         % seconds) divided by the number of slices, timing(1) = time between slices
%         timings(2) = timings(1)+(1000/1000); % time to acquire the last slice - our sequence
%         % includes a gap between the end of the last slice and the beginning of the next volume
%         % timing(2) = time between last slice and next volume
%
%
%         for e=1:size(my_epi,2)
%         % do the slice timing one session at a time - we don't want to
%         % interpolate across scans from different sessions
%             files = get_files(my_epi{e},'rf*.nii');
%             spm_slice_timing(files, sliceorder, refslice, timings);
%         end
%
%
%         % this step should create a new set of images with the prefix 'arf*.nii'
%     else
%         disp('skipping step 4 slice timing as output file exists');
%         % files is out of date - but not app used again
%     end

% =========================================================================


% =========================================================================
% Step  5 - Coregister the structural to the mean undistorted EPI
% =========================================================================

coreg = 1;

if coreg
    
    disp('running co-reg struct');
    
    my_struc = fullfile(rootdir, 'Results_temp_playing', csub, struc);
    mean_img = get_files(my_epi{1},'meanf*.nii'); % mean undistorted image
    struc_img = get_files(my_struc,'s*.nii'); % trimmed structural image
    
    options.cost_fun = 'nmi'; % which cost function to use - normalised mutual information
    options.sep = [4 2]; % resolution at which to sample images
    options.tol = [0.0200    0.0200    0.0200    0.0010    0.0010    0.0010 ...
        0.0100    0.0100 0.0100    0.0010    0.0010    0.0010]; % tolerance for parameters
    options.fwhm = [7 7]; % degree of smoothing to apply to images before registration
    
    x=spm_coreg(mean_img,struc_img,options); % find the coreg parameters
    
    % convert the coreg parameters into an affine transformation matrix
    M  = inv(spm_matrix(x));
    MM = zeros(4,4,1);
    MM(:,:,1) = spm_get_space(struc_img);
    
    % modify the header of the structural image
    spm_get_space(struc_img, M*MM(:,:,1));
    
    % print the results screen
    %spm_print(printfile);
    spm_print(fullfile(rootdir,'preprocessing',['spm_coregister_',num2str(csub,'S%03d'),'.ps']));
    
else
    disp('skipping coregistration');
end



% =========================================================================


% =========================================================================
% Step  6 - Normalise the structural
% =========================================================================
%
% % The original method of doing this was to manually trim and skull strip
% the structural (using trim_img and bet), then to match the template using
% the whole image. The new method, which we're going to use here, is more
% automatic, and involves segmenting the image into white and grey matter
% first, then normalising these to tissue specific templates. This avoids
% the need for skull stripping.

% It does however benefit from a 2 pass procedure... The first pass doesn't
% try to segment the image, it just corrects any bias present in the image
% - i.e. any systematic difference in signal between different parts of the
% brain. This is necessary as the MPRAGE images tend to be darker at the
% front than the back

normalise_struct = 1;

if normalise_struct
    
    disp('running normalise struc');
    struc_img = get_files(my_struc,'s*.nii'); % trimmed structural image
    
    
    %%%%%%%% 1st pass:
    % -------------------------------------------------------------------------
    estopts.regtype='';                     % turn off affine registration
    out = spm_preproc(struc_img,estopts);   % estimate bias field
    sn = spm_prep2sn(out);                  % convert to a transformation matrix
    
    writeopts.biascor = 1;                  % only write out attenuation corrected image
    writeopts.GM  = [0 0 0];                % turn off everything else...
    writeopts.WM  = [0 0 0];
    writeopts.CSF = [0 0 0];
    writeopts.cleanup = 0;
    spm_preproc_write(sn,writeopts);        % write bias corrected image (prepends 'm' suffix)
    
    
    
    %%%%%%%% 2nd pass using attenuation corrected image
    % -------------------------------------------------------------------------
    struc_img = get_files(my_struc,'ms*.nii'); % corrected structural image
    
    estopts.regtype='mni';    % turn on affine again
    out = spm_preproc(struc_img,estopts);       % estimate normalisation parameters
    [sn,isn] = spm_prep2sn(out);                % convert to matrix
    
    % write out GM and WM native + unmod normalised
    writeopts.biascor = 1;
    writeopts.GM  = [0 1 1];
    writeopts.WM  = [0 1 1];
    writeopts.CSF = [0 0 0];
    writeopts.cleanup = 0;
    spm_preproc_write(sn,writeopts);
    
    % save normalisation parametrs to a matrix file - these will be used at the
    % next stage to normalise the functional data. An inverse matrix file is
    % also created, which is useful if you want to "un-normalise" anything  -
    % e.g. regions of interest that are defined in template space.
    [pth fle]=fileparts(struc_img);
    matname = fullfile(pth,[fle '_seg_sn.mat']);
    invmatname = fullfile(pth,[fle '_seg_inv_sn.mat']);
    %savefields(matname,sn);
    %savefields(invmatname,isn);
    
    % Save fields
    fn = fieldnames(sn);
    for i=1:length(fn)
        eval([fn{i} '= sn.' fn{i} ';']);
    end
    if str2double(version('-release'))>=14
        save(matname,'-V6',fn{:});
    else
        save(matname,fn{:});
    end
    fn_inv = fieldnames(isn);
    for i=1:length(fn)
        eval([fn{i} '= isn.' fn{i} ';']);
    end
    if str2double(version('-release'))>=14
        save(invmatname,'-V6',fn_inv{:});
    else
        save(invmatname,fn_inv{:});
    end
    
    spm_write_sn(struc_img,matname); % write out normalised structural - this
    % is only really to check that the normalisation has worked, its not used
    % in any further analyses, so its only written out at fairly low resolution.
    
    
    %%%%%%%% Display the results
    % -------------------------------------------------------------------------
    try
        figure(spm_figure('FindWin'));
    end
    def.temp = fullfile([Main_analysis_directory,'\spm12\canonical\avg152T1.nii']);

%     def.temp = 'F:/Postdoc/MD_project/fMRI/Data/Analyses/spm12/canonical/avg152T1.nii';
    
    cd(fullfile(rootdir,'Results_temp_playing',num2str(csub,'S%03d'),struc));
    imgs = char(def.temp,...  % T1 template
        struc_img,...                                           % Un-normalised structural
        char(get_files(my_struc,'wms*.nii')),...                      % Normalised structural
        char(get_files(my_struc,'c1m*.nii')));                        % Un-normalised grey matter
    
    imgs = spm_vol(imgs);
    spm_check_registration(imgs); % display the images
    
    ann1=annotation('textbox',[.1 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','T1 template');
    ann2=annotation('textbox',[.6 .891 .3 .025],'HorizontalAlignment','center','Color','r','String','Native T1');
    ann3=annotation('textbox',[.1 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Normalised T1');
    ann4=annotation('textbox',[.6 .413 .3 .025],'HorizontalAlignment','center','Color','r','String','Native segmented grey matter');
    
    %spm_print(printfile);
    spm_print(fullfile(rootdir,'preprocessing',['spm_normalisation_',num2str(csub,'S%03d'),'.ps']));
    
    spm_print(pdffile);
    
    delete(ann1); delete(ann2); delete(ann3); delete(ann4); f=figure(spm_figure('FindWin')); clf(f);
    
else
    disp('skipping normalise structural');
end

% =========================================================================


%     % =========================================================================
%     % Step  7 - Apply the normalisation parameters to the EPI images
%     % =========================================================================
%

normaliseEPIS = 1; %turn on/off normalise EPIS here

if normaliseEPIS
    % check if warped files exist
    %         wfiles = [];
    %         wfiles = get_files(char(my_epi),'warf*.nii');
    %
    %         if isempty(wfiles) %normalise if don't exist
    
    disp('running normalise EPI');
    options.bb = [  -78  -112   -50; 78    76    85];   % = bounding box = the range of
    % co-ordinates (in mm) to include in the image
    options.vox = [2 2 2];                              % size of voxels to use in the normalised images
    options.interp = 1;                                 % interpolation method
    options.wrap = [0 0 0];                             % wrap edges?
    options.preserve = 0;                               % preserve voxel concentrations? (mainly for VBM -
    % see spm_write_sn for more details)
    
    files = get_files(char(my_epi),'arf*.nii'); % Undistorted EPIs
    files = char(files, get_files(my_epi{1},'meanf*.nii')); % Add mean undistorted EPI
    
    tic();
    spm_write_sn(files,matname,options);
    toc();
    
else
    disp('skipping normalise EPI');
end



% =========================================================================
% Step  8 - Smooth the normalised EPI images
% =========================================================================

smoothnorm = 1;

if smoothnorm == 1;
    
    disp('running smooth norm epis');
    tic();
    
    files = get_files(char(my_epi),'warf*.nii'); % EPI images
    files = char(files,get_files(my_epi{1},'wmeanf*.nii')); % add mean undistorted EPI
    
    %         smoothing_name = 's8';
    smoothing_name = 's4';
    
    % check if smoothed files exist
    sfiles = [];
    sfiles = get_files(char(my_epi),'s4*.nii');
    
    %if isempty(sfiles) % only smooth if files don't exist
    
    %             FWHM = [8 8 8]; % amount of smoothing to apply - fwhm of smoothing kernel in [x,y,z] in mm
    FWHM = [4 4 4]; % amount of smoothing to apply - fwhm of smoothing kernel in [x,y,z] in mm
    
    
    % Have to loop this manually as spm_smooth only does one image at a time...
    n = size(files,1);
    spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
    
    for i = 1:n
        cimg = deblank(files(i,:)); % current image
        [pth,fle,ext] = fileparts(cimg);
        simg = fullfile(pth,[smoothing_name fle ext]);
        %simg = fullfile(pth,['s4' fle ext]); % create output file name - prepend 's10' suffix
        %if ~exist(simg)
        spm_smooth(cimg,simg,FWHM); % do the smooth
        %end
        spm_progress_bar('Set',i); % update the progress bar
    end
    spm_progress_bar('Clear');
    
    %disp(['finished ' csub]);
    toc();
    
    %end
else
    disp('skipping smooth normalised');
end

% =========================================================================
% Step  9 - Smooth the non normalised EPI images
% =========================================================================


smoothnative = 1;

if smoothnative == 1
    
    disp('running smooth non norm epis');
    tic();
    
    files = get_files(char(my_epi),'arf*.nii'); % EPI images
    files = char(files,get_files(my_epi{1},'meanf*.nii')); % add mean undistorted EPI
    
    smoothing_name = 's4';
    
    FWHM = [4 4 4]; % amount of smoothing to apply - fwhm of smoothing kernel in [x,y,z] in mm
    
    
    % Have to loop this manually as spm_smooth only does one image at a time...
    n = size(files,1);
    spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
    
    for i = 1:n
        cimg = deblank(files(i,:)); % current image
        [pth,fle,ext] = fileparts(cimg);
        simg = fullfile(pth,[smoothing_name fle ext]);
        %simg = fullfile(pth,['s4' fle ext]); % create output file name - prepend 's10' suffix
        %if ~exist(simg)
        spm_smooth(cimg,simg,FWHM); % do the smooth
        %end
        spm_progress_bar('Set',i); % update the progress bar
    end
    spm_progress_bar('Clear');
    
    disp(['finished ' csub]);
    toc();
    
    %     catch
    %         disp(['processing failed:  ' csub]);
    %         failedsubs = [failedsubs csub];
    %     end
    
else
    disp('skipping smooth non normalised');
end

%print warnings
% if failedsubs;
%     disp(['THESE SUBJECTS FAILED TO PROCESS: ' failedsubs]);
% end
% =========================================================================

% convert ps to pdf
%ps2pdf('psfile', printfile, 'pdffile', pdffile);


% =========================================================================
% =========================================================================
% Subfunctions
% =========================================================================
% =========================================================================






% =========================================================================
    function savefields(fnam,p)
        % =========================================================================
        
        if length(p)>1, error('Can''t save fields.'); end;
        fn = fieldnames(p);
        if numel(fn)==0, return; end;
        for i=1:length(fn),
            eval([fn{i} '= p.' fn{i} ';']);
        end;
        if str2double(version('-release'))>=14,
            save(fnam,'-V6',fn{:});
        else
            save(fnam,fn{:});
        end;
        
        % return;
    end
done=1;
end