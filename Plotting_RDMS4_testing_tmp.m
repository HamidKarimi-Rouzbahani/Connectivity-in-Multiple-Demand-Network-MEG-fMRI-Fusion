clc;
close all;
clear all;
Main_analysis_directory='D:/Hamid/Postdoc/MD/Analyses';

Action=4;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs

Trial_type=1; % 1=correct (for RSA=colors*rules*stimuli); 2=same errors (for RSA= the same as previous);
% 3=swapped errors (only for decoding)
Stim_Resp='Stim';

% 1 (perfect)= columns: rule1inner rule1outer rule2inner rule2outer
% 5=peripheral stimulus coding (only for RSA: collapsed across stimulus side and cue colors from different rules)

% 4=stimulus coding  (only for RSA: 1 to 4 stimuli ):

% 6=side stimulus coding (the same as 5; only for RSA: collapsed across cue colors from different rules)
% 10 (perfect)=Stim side coding (with rules included: columns are rule1stim1 rule1stim2 rule2stim1 rule2stim2)

% 3 (perfect)=rule coding across cue colors (only for RSA: 1 to 4 cue colors collapsed across all stimuli)

% 7=Response coding (responses 1-4): bad because stimuli overwhelmed responses
% 8=Response coding (with rules included: columns are
% rule1resp1(leftfingers) rule1resp2(rightfingers) rule2resp1 rule2resp2): not good as it makes
% nagative corrs in frontal areas where rules are coded
% 9=Response coding (with rules included: columns are color1resp1(leftfingers) color1resp2(rightfingers)
% color2resp1 color2resp2)
% 11 (perfect)=Response coding Inner vs Outer and not left and right(with rules included: columns are color1resp1 color1resp2 color2resp1 color2resp2)
% 12 (perfect)=same as 11 but different colours collapsed


smoothing_rate=100;


Decod_instead_of_Corr=1;
searchlighting=0;

subjects=[1];

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


addpath('D:\Hamid\Postdoc\MD\MEG')
if Trial_type==1
    %
    %     load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4.mat']);
    %     sampling_time=25;
    
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4_5ms.mat']);
    sampling_time=5;
    
    
    
    
elseif Trial_type==3
    %     load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_cue_coding_RDM.mat']);
    %     sampling_time=5;
    
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_cue_coding_RDM_dec.mat']);
    sampling_time=5;
    
elseif Trial_type==4
    %     load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_stm_coding_RDM.mat']);
    %     sampling_time=5;
    
    % the following is broken at the moment
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_stm_coding_RDM_dec.mat']);
    sampling_time=25;
    
elseif Trial_type==5
    
    %     load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Stm_coding_Clpsd_cue_colors_RDM.mat']);
    %     sampling_time=5;
    
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Stm_coding_Clpsd_cue_colors_RDM_dec.mat']);
    sampling_time=25;
    
elseif Trial_type==6
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Stm_side_coding_Clpsd_cue_colors_RDM_dec.mat']);
    sampling_time=25;
    sampling_time=5;
    
    
elseif Trial_type==7
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==8
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_acrs_rules_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==9
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_acrs_opst_col_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==10
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4_5ms_side_dec.mat']);
    sampling_time=5;
    
elseif Trial_type==11
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_InOut_acrs_opst_col_RDM_decoding.mat']);
    sampling_time=5;
elseif Trial_type==12
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_InOut_acrs_opst_col_RDM_decoding.mat']);
    sampling_time=5;
end

for i=1:4
    RDM_correct_4_all(i,i,:)=nan;
    RDM_error_4_all(i,i,:)=nan;
end

Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
MEG_RDM=RDM_correct_4_all;
NoC=4;
if Trial_type==3
    Beta_folder_name=[Beta_folder_name,'_CueCoding'];
elseif Trial_type==4
    Beta_folder_name=[Beta_folder_name,'_StmCoding'];
elseif Trial_type==5
    Beta_folder_name=[Beta_folder_name,'_StmPeriCod_clpsd_OppstCuCol'];
elseif Trial_type==6
    Beta_folder_name=[Beta_folder_name,'_StmSideCod_clpsd_OppstCuCol'];
elseif Trial_type==7
    Beta_folder_name=[Beta_folder_name,'_RespCod_1to4'];
elseif Trial_type==8
    Beta_folder_name=[Beta_folder_name,'_RespCod_inRules'];
elseif Trial_type==9
    Beta_folder_name=[Beta_folder_name,'_RespCod_clpsd_OppstCuCol'];
elseif Trial_type==10
    Beta_folder_name=[Beta_folder_name,'_StimSideCoding'];
elseif Trial_type==11
    Beta_folder_name=[Beta_folder_name,'_RespCod_InOut_clpsd_OppstCuCol'];
elseif Trial_type==12
    Beta_folder_name=[Beta_folder_name,'_RespCod_InOut_clpsd2_OppstCuCol'];
end

% Regions={'ACC SMA','Left DLPFC','Left IPS','Left VLPFC',...
%     'Right DLPFC','Right IPS','Right VDPFC'...
%     'Left LOC','Right LOC','Left Visual','Right Visual'};
regions_order_to_read=[1 3 6 2 5 4 7 8 9 10 11];

Regions={'ACC SMA','Left IPS','Right IPS','Left DLPFC','Right DLPFC',...
    'Left VLPFC','Right VLPFC','Left LOC','Right LOC','Left Visual','Right Visual'};

Regions_summereized={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};
Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];


for region=[1:11]
    
    c=0;
    for subject=subjects
        c=c+1;
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        tmp_RSA(:,:,c)=nan*ones(NoC);
        if Trial_type==1
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Dec_ROIs.mat']),'RSA_results');
        elseif Trial_type==3
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Cue_Color_Dec_ROIs.mat']),'RSA_results');
        elseif Trial_type==4
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stm_Dec_ROIs.mat']),'RSA_results');
        elseif Trial_type==5
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stm_clpsdCueCol_Dec_ROIs.mat']),'RSA_results');
        elseif Trial_type==6
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stm_Side_clpsdCueCol_Dec_ROIs.mat']),'RSA_results');
        elseif Trial_type==7
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_Dec_ROIs.mat']),'RSA_results');
        elseif Trial_type==8
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_acros_CuesDec_ROIs.mat']),'RSA_results');
        elseif Trial_type==9
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_opst_ColrsDec_ROIs.mat']),'RSA_results');
        elseif Trial_type==10
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stim_side_acrs_rulesDec_ROIs.mat']),'RSA_results');
        elseif Trial_type==11
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_InOut_opst_ColrsDec_ROIs.mat']),'RSA_results');
        elseif Trial_type==12
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_InOut_opst_Colrs2Dec_ROIs.mat']),'RSA_results');
        end
        combinations=nchoosek([1:NoC],2);
        for i=1:length(combinations)
            if isempty(RSA_results{region, i})
                tmp_RSA(combinations(i,1),combinations(i,2),c)=nan;
            else
                tmp_RSA(combinations(i,1),combinations(i,2),c)=RSA_results{region, i,1}.accuracy.output;
            end
        end
        
        for i=1:size(tmp_RSA,1)
            for j=1:size(tmp_RSA,2)
                if i>=j
                    tmp_RSA(i,j,c)=nan;
                end
            end
        end
        
        for i=1:length(combinations)
            g=0;
            for iter=2:5
                g=g+1;
                tmp_RSA_rand(combinations(i,1),combinations(i,2),g)=RSA_results{region,i,iter}.accuracy.output;
            end
        end
        
        for i=1:size(tmp_RSA,1)
            for j=1:size(tmp_RSA,2)
                if i>=j
                    tmp_RSA_rand(i,j,:)=nan;
                end
            end
        end
        
    end        
end
subplot(121)
imagesc(tmp_RSA,[50 100])
subplot(122)