clc;
close all;
clear all;
Main_analysis_directory='Z:/projects/Hamid/Projects/MD/Analyses';

Action=4;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs
tt=0;
Stim_Resp='Stim';
datas={'StimSide','StimPeriph','Rule','Response'};
% for Trial_type=[10 1 3 11]
Trial_type=[11]
% 3=swapped errors (only for decoding)

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




Decod_instead_of_Corr=1;
searchlighting=0;

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



addpath('Z:/projects/Hamid/Projects/MD/MEG')
if Trial_type==1
    %
    %     load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4.mat']);
    %     sampling_time=25;
    
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4_5ms.mat']);
    sampling_time=5;
    
    
    
    
elseif Trial_type==3
    %     load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_cue_coding_RDM.mat']);
    %     sampling_time=5;
    
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_cue_coding_RDM_dec.mat']);
    sampling_time=5;
    
elseif Trial_type==4
    %     load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_stm_coding_RDM.mat']);
    %     sampling_time=5;
    
    % the following is broken at the moment
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_stm_coding_RDM_dec.mat']);
    sampling_time=25;
    
elseif Trial_type==5
    
    %     load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Stm_coding_Clpsd_cue_colors_RDM.mat']);
    %     sampling_time=5;
    
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Stm_coding_Clpsd_cue_colors_RDM_dec.mat']);
    sampling_time=25;
    
elseif Trial_type==6
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Stm_side_coding_Clpsd_cue_colors_RDM_dec.mat']);
    sampling_time=25;
    sampling_time=5;
    
    
elseif Trial_type==7
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==8
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_acrs_rules_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==9
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_acrs_opst_col_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==10
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4_5ms_side_dec.mat']);
    sampling_time=5;
    
elseif Trial_type==11
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_InOut_acrs_opst_col_RDM_decoding.mat']);
    sampling_time=5;
    
elseif Trial_type==12
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_InOut_acrs_opst_col2_RDM_decoding.mat']);
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


Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];

for region=[1:11]
    
    c=0;
    for subject=subjects
        c=c+1;
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        tmp_RSA(:,:,c)=nan(NoC);
%         if Trial_type==1
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Dec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==3
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Cue_Color_Dec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==4
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stm_Dec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==5
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stm_clpsdCueCol_Dec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==6
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stm_Side_clpsdCueCol_Dec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==7
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_Dec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==8
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_acros_CuesDec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==9
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_opst_ColrsDec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==10
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Stim_side_acrs_rulesDec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==11
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_InOut_opst_ColrsDec_ROIs.mat']),'RSA_results');
%         elseif Trial_type==12
%             load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Iter_Resp_InOut_opst_Colrs2Dec_ROIs.mat']),'RSA_results');
%         end
            if Trial_type==1
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Dec_ROIs.mat']),'RSA_results');
            elseif Trial_type==3
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Cue_Color_Dec_ROIs.mat']),'RSA_results');
            elseif Trial_type==4
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Stm_Dec_ROIs.mat']),'RSA_results');
            elseif Trial_type==5
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Stm_clpsdCueCol_Dec_ROIs.mat']),'RSA_results');
            elseif Trial_type==6
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Stm_Side_clpsdCueCol_Dec_ROIs.mat']),'RSA_results');
            elseif Trial_type==7
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Resp_Dec_ROIs.mat']),'RSA_results');
            elseif Trial_type==8
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Resp_acros_CuesDec_ROIs.mat']),'RSA_results');
            elseif Trial_type==9
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Resp_opst_ColrsDec_ROIs.mat']),'RSA_results');
            elseif Trial_type==10
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Stim_side_acrs_rulesDec_ROIs.mat']),'RSA_results');
            elseif Trial_type==11
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Resp_InOut_opst_ColrsDec_ROIs.mat']),'RSA_results');
            elseif Trial_type==12
                load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Resp_InOut_opst_Colrs2Dec_ROIs.mat']),'RSA_results');
             end
        combinations=nchoosek([1:NoC],2);
        for i=1:length(combinations)
            if isempty(RSA_results{region, i,1})
                tmp_RSA(combinations(i,1),combinations(i,2),c)=nan;
            else
                tmp_RSA(combinations(i,1),combinations(i,2),c)=RSA_results{region,i,1}.accuracy.output;
            end
        end
        for i=1:size(tmp_RSA,1)
            for j=1:size(tmp_RSA,2)
                if i>=j
                    tmp_RSA(i,j,c)=nan;
                end
            end
        end
    end
    RSA_fMRI(:,:,region)=nanmean(tmp_RSA,3);
end
load(['Z:/projects/Hamid/Projects/MD/MEG/RDM_models_',num2str(NoC),'.mat'])


%% Plotting
% fMRI-model correlation
smoothing_rate=20;
% Regions={'ACC SMA','Left DLPFC','Left IPS','Left VLPFC',...
%     'Right DLPFC','Right IPS','Right VLPFC'...
%     'Left LOC','Right LOC','Left Visual','Right Visual'};
regions_order_to_read=[11 10 9 8 6 3 7 4 5 2 1];
Regions_summereized={'Visual','LOC','IPS','AI/FO','IFS','ACC'};

RDM_fMRIs(:,1)=mean([reshape(RSA_fMRI(:,:,11),[NoC*NoC 1]) reshape(RSA_fMRI(:,:,10),[NoC*NoC 1])],2);
RDM_fMRIs(:,2)=mean([reshape(RSA_fMRI(:,:,9),[NoC*NoC 1]) reshape(RSA_fMRI(:,:,8),[NoC*NoC 1])],2);
RDM_fMRIs(:,3)=mean([reshape(RSA_fMRI(:,:,6),[NoC*NoC 1]) reshape(RSA_fMRI(:,:,3),[NoC*NoC 1])],2);
RDM_fMRIs(:,4)=mean([reshape(RSA_fMRI(:,:,7),[NoC*NoC 1]) reshape(RSA_fMRI(:,:,4),[NoC*NoC 1])],2);
RDM_fMRIs(:,5)=mean([reshape(RSA_fMRI(:,:,5),[NoC*NoC 1]) reshape(RSA_fMRI(:,:,2),[NoC*NoC 1])],2);
RDM_fMRIs(:,6)=reshape(RSA_fMRI(:,:,1),[NoC*NoC 1]);

% RDM_fMRIs(:,1)=mean([reshape(RSA_fMRI(:,:,11),[NoC*NoC 1])],2);
% RDM_fMRIs(:,2)=mean([reshape(RSA_fMRI(:,:,9),[NoC*NoC 1])],2);
% RDM_fMRIs(:,3)=mean([reshape(RSA_fMRI(:,:,6),[NoC*NoC 1])],2);
% RDM_fMRIs(:,4)=mean([reshape(RSA_fMRI(:,:,7),[NoC*NoC 1])],2);
% RDM_fMRIs(:,5)=mean([reshape(RSA_fMRI(:,:,5),[NoC*NoC 1])],2);
% RDM_fMRIs(:,6)=reshape(RSA_fMRI(:,:,1),[NoC*NoC 1]);

figure;
subplot(221)
if Trial_type==10
    titles={'StimSide'};
elseif Trial_type==1
    titles={'StimPeriph'};
elseif Trial_type==3
    titles={'Rule'};
elseif Trial_type==11
    titles={'Response'};
end
c=0;
for region=1:6
    c=c+1;
    RDM_fMRI=RDM_fMRIs(~isnan(RDM_fMRIs(:,region)),region);
    if Trial_type==1 || Trial_type>3
        RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
        RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
        cors(c)=corr(RDM_fMRI,RDM_model_Stim,'type','Pearson');
    elseif Trial_type==3
        RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
        RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));
        cors(c)=corr(RDM_fMRI,RDM_model_Rule,'type','Pearson');
    end
    
end

bar([1:size(cors,2)],cors)
hold on
xticks([1:length(Regions_summereized)])
xticklabels(Regions_summereized)
xtickangle(45)
ylabel('fMRI-model correlation')
title(titles)
box off

% figure;
subplot(223)
cors_inter_area=nan(6);
for region=1:6
    RDM_fMRI_r1=RDM_fMRIs(~isnan(RDM_fMRIs(:,region)),region);
    for region2=region+1:6
        RDM_fMRI_r2=RDM_fMRIs(~isnan(RDM_fMRIs(:,region2)),region2);
        [cors_inter_area(region2,region),ps(region2,region)]=corr(RDM_fMRI_r1,RDM_fMRI_r2,'type','Pearson');
    end
end
imagesc(cors_inter_area)
hold on;
% correct for MC
% ps(~isnan(ps))=mafdr(ps(~isnan(ps)));
for i=1:6
    for j=i+1:6
        if ps(j,i)<0.05
            text(i-.3,j,[num2str(cors_inter_area(j,i),'%1.2f'),'*'],'FontSize',15);
        else
            text(i-.3,j,num2str(cors_inter_area(j,i),'%1.2f'),'FontSize',15);
        end
    end
end
if Trial_type==10
    titles={'RoI-RoI correlation (StimSide)'};
elseif Trial_type==1
    titles={'RoI-RoI correlation (StimPeriph)'};
elseif Trial_type==3
    titles={'RoI-RoI correlation (Rule)'};
elseif Trial_type==11
    titles={'RoI-RoI correlation (Response)'};
end
xticks([1:6])
xticklabels(Regions_summereized)
yticks([1:6])
yticklabels(Regions_summereized)
title(titles)
xtickangle(45)
ytickangle(45)
box off




% MEG-model correlation
% figure;
subplot(222)
MEG_RDM=RDM_correct_4_all;

if strcmp(Stim_Resp,'Stim')
    Times=[-500:sampling_time:3995]+sampling_time./2;
    xlim([-200 2500])
else
    Times=[-4000:sampling_time:995]+sampling_time./2;
    xlim([-2500 500])
end

for time=1:length(Times)
    if Trial_type==1 || Trial_type>3
        correl_MEG(time)=corr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(S_model_cor,[NoC*NoC 1]),'type','Pearson','rows','complete');
    elseif Trial_type==3
        correl_MEG(time)=corr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(R_model_cor,[NoC*NoC 1]),'type','Pearson','rows','complete');
    end
end
plots=plot(Times,smooth(correl_MEG,smoothing_rate),'linewidth',2);

if strcmp(Stim_Resp,'Stim')
    xlim([-200 2500])
else
    xlim([-2500 500])
end
line([min(Times) max(Times)],[0 0],'linestyle','--','color','k')
line([0 0],[-0.6 0.8],'linestyle','--','color','k')
xlabel('Time (ms)')
ylabel('MEG-model correlation')
box off
if Trial_type==10
    titles={'StimSide'};
elseif Trial_type==1
    titles={'StimPeriph'};
elseif Trial_type==3
    titles={'Rule'};
elseif Trial_type==11
    titles={'Response'};
end
title (titles)



% fMRI-MEG correlation
% figure;
subplot(224)
for region=1:6
    for time=1:length(Times)
        RDM_MEG=reshape(MEG_RDM(:,:,time),[NoC*NoC 1]);
        RDM_MEG=RDM_MEG(~isnan(RDM_MEG));
        RDM_fMRI=RDM_fMRIs(:,region);
        RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
        MEG_fMRI(time,region)=corr(RDM_MEG,RDM_fMRI,'type','Pearson');
        
        [Trial_type region time]
    end
end

for region=1:6
    plot(Times,smooth(MEG_fMRI(:,region),smoothing_rate),'linewidth',2)
    hold on;
end

if strcmp(Stim_Resp,'Stim')
    xlim([-200 2500])
else
    xlim([-2500 500])
end
line([min(Times) max(Times)],[0 0],'linestyle','--','color','k')
line([0 0],[-0.6 0.8],'linestyle','--','color','k')
xlabel('Time (ms)')
ylabel('MEG-fMRI correlation')
box off
if Trial_type==10
    titles={'StimSide'};
elseif Trial_type==1
    titles={'StimPeriph'};
elseif Trial_type==3
    titles={'Rule'};
elseif Trial_type==11
    titles={'Response'};
end
title (titles)
legend(Regions_summereized)
