clc;
close all;
clear all;
Main_analysis_directory='Z:/projects/Hamid/Projects/MD/Analyses';

Action=3;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs
Trial_type=1; % 1=correct; 2=all errors

Stim_Resp='Resp';

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
if Action==3
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding.mat']);
elseif Action==4
    load(['Z:/projects/Hamid/Projects/MD/MEG/Subjects_Average_',Stim_Resp,'_aligned_RDM_decoding4.mat']);
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
        load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Dec_ROIs.mat']),'RSA_results');
        combinations=nchoosek([1:NoC],2);
        for i=1:length(combinations)
            if isempty(RSA_results{region, i})
                tmp_RSA(combinations(i,1),combinations(i,2),c)=nan;
            else
                tmp_RSA(combinations(i,1),combinations(i,2),c)=RSA_results{region, i}.accuracy.output;
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
load(['Z:/projects/Hamid/Projects/MD/MEG\RDM_models_',num2str(NoC),'.mat'])
%% Plotting fMRI and mean MEG RDMs
scales=[nanmin(nanmin(nanmin(RSA_fMRI))) nanmax(nanmax(nanmax(RSA_fMRI)))];
c=0;
for region=regions_order_to_read
    c=c+1;
    subplot_tmp=subplot(4,4,c);
    hold(subplot_tmp,'on');
    image(RSA_fMRI(:,:,region),'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    
    xticks(1:NoC)
    yticks(1:NoC)
    ytickangle(45)
    xtickangle(45)
    title(Regions{c})
end
subplot_tmp=subplot(4,4,12);
hold(subplot_tmp,'on');
if Action==2
    image(nanmean(RDM_correct_16_all,3),'Parent',subplot_tmp,'CDataMapping','scaled');
    scales=[nanmin(nanmin(nanmean(RDM_correct_16_all,3))) nanmax(nanmax(nanmean(RDM_correct_16_all,3)))];
elseif Action==3
    image(nanmean(RDM_correct_8_all,3),'Parent',subplot_tmp,'CDataMapping','scaled');
    scales=[nanmin(nanmin(nanmean(RDM_correct_8_all,3))) nanmax(nanmax(nanmean(RDM_correct_8_all,3)))];
elseif Action==4
    image(nanmean(RDM_correct_4_all,3),'Parent',subplot_tmp,'CDataMapping','scaled');
    scales=[nanmin(nanmin(nanmean(RDM_correct_4_all,3))) nanmax(nanmax(nanmean(RDM_correct_4_all,3)))];
end
axis(subplot_tmp,'tight');
axis(subplot_tmp,'ij');
set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
xticks(1:NoC)
yticks(1:NoC)
ytickangle(45)
xtickangle(45)
title('Mean of MEG')

subplot_tmp=subplot(4,4,13);
hold(subplot_tmp,'on');
image(S_model_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
axis(subplot_tmp,'tight');
axis(subplot_tmp,'ij');
set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
title('Stimulus Model')

subplot_tmp=subplot(4,4,14);
hold(subplot_tmp,'on');
image(R_model_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
axis(subplot_tmp,'tight');
axis(subplot_tmp,'ij');
set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
title('Rule Model')

subplot_tmp=subplot(4,4,15);
hold(subplot_tmp,'on');
image(S_model_periph_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
axis(subplot_tmp,'tight');
axis(subplot_tmp,'ij');
set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
title('Stim Side')

subplot_tmp=subplot(4,4,16);
hold(subplot_tmp,'on');
image(S_model_decod_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
axis(subplot_tmp,'tight');
axis(subplot_tmp,'ij');
set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
title('Stim Periph')

%% plotting fMRI correlation to model RDMs
figure;
titles={'Stimulus side','Stimulus','Stimulus periph','Rule'};
c=0;
for region=regions_order_to_read
    c=c+1;
    RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
    RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
    
    RDM_model_Stim_side=reshape(S_model_periph_cor,[NoC*NoC 1]);
    RDM_model_Stim_side=RDM_model_Stim_side(~isnan(RDM_model_Stim_side));
    
    RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
    RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
    
    RDM_model_Stim_periph=reshape(S_model_decod_cor,[NoC*NoC 1]);
    RDM_model_Stim_periph=RDM_model_Stim_periph(~isnan(RDM_model_Stim_periph));

    RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
    RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));

    cors(1,c)=partialcorr(RDM_fMRI,RDM_model_Stim_side,[RDM_model_Stim RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');    
    cors(2,c)=partialcorr(RDM_fMRI,RDM_model_Stim,[RDM_model_Stim_side RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');
    cors(3,c)=partialcorr(RDM_fMRI,RDM_model_Stim_periph,[RDM_model_Stim_side RDM_model_Stim RDM_model_Rule],'type','Spearman');
    cors(4,c)=partialcorr(RDM_fMRI,RDM_model_Rule,[RDM_model_Stim RDM_model_Stim_side RDM_model_Stim_periph],'type','Spearman');
end
cors_summerized=[cors(:,1) nanmean(cors(:,2:3),2) nanmean(cors(:,4:5),2),...
    nanmean(cors(:,6:7),2) nanmean(cors(:,8:9),2) nanmean(cors(:,10:11),2)];
for info=1:4
%     subplot(2,2,info);
%     bar([1:size(cors,2)],cors(info,:))
%     hold on
%     xticks([1:length(Regions)])
%     xticklabels(Regions)
%     xtickangle(45)
%     ylim([-0.4 0.8])
%     ylabel('Spearman''s correlation to model(\rho)')
%     title(titles{info})
%     box off
%     
    subplot(2,2,info);
    bar([1:size(cors_summerized,2)],cors_summerized(info,:))
    hold on
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
    xtickangle(45)
    ylim([-0.4 0.8])
    ylabel('Correlation (\rho)')
    title(titles{info})
    box off
end
ccc
%% Plotting MEG correlation to model RDMs
% close all;
figure;
if NoC==16
    MEG_RDM=RDM_correct_16_all;
    MEG_RDM_err=RDM_error_4_all;
elseif NoC==8
    MEG_RDM=RDM_correct_8_all;
    MEG_RDM_err=RDM_error_8_all;
elseif NoC==4
    MEG_RDM=RDM_correct_4_all;
    MEG_RDM_err=RDM_error_16_all;
end

if strcmp(Stim_Resp,'Stim')
    Times=[-500:25:3995];
else
    Times=[-4000:25:995];
end

for time=1:length(Times)
    
    RDM_model_Stim_side=reshape(S_model_periph_cor,[NoC*NoC 1]);
    RDM_model_Stim_side=RDM_model_Stim_side(~isnan(RDM_model_Stim_side));
    
    RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
    RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
    
    RDM_model_Stim_periph=reshape(S_model_decod_cor,[NoC*NoC 1]);
    RDM_model_Stim_periph=RDM_model_Stim_periph(~isnan(RDM_model_Stim_periph));

    MEG_rdm=reshape(MEG_RDM(:,:,time),[NoC*NoC 1]);
    MEG_rdm=MEG_rdm(~isnan(MEG_rdm));
    
    correl(1,time)=partialcorr(MEG_rdm,RDM_model_Stim_side,[RDM_model_Stim RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');    
    correl(2,time)=partialcorr(MEG_rdm,RDM_model_Stim,[RDM_model_Stim_side RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');
    correl(3,time)=partialcorr(MEG_rdm,RDM_model_Stim_periph,[RDM_model_Stim_side RDM_model_Stim RDM_model_Rule],'type','Spearman');
    correl(4,time)=partialcorr(MEG_rdm,RDM_model_Rule,[RDM_model_Stim RDM_model_Stim_side RDM_model_Stim_periph],'type','Spearman');
    correl(5,time)=corr(MEG_rdm,RDM_model_Rule,'type','Spearman');
    
%     correl(1,time)=partialcorr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(S_model_periph_cor,[NoC*NoC 1]),[reshape(R_model_cor,[NoC*NoC 1]) reshape(R_model_cor,[NoC*NoC 1]) reshape(R_model_cor,[NoC*NoC 1])],'type','Spearman','rows','complete');
%     correl(2,time)=partialcorr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(R_model_cor,[NoC*NoC 1]),[reshape(R_model_cor,[NoC*NoC 1]) reshape(R_model_cor,[NoC*NoC 1]) reshape(S_model_periph_cor,[NoC*NoC 1])],'type','Spearman','rows','complete');

end


for i=1:5
    plots(i)=plot(Times,smooth(correl(i,:),10),'linewidth',2);
    hold on;
end
legend ([plots(1) plots(2) plot(3) plot(4)],titles)
line([min(Times) max(Times)],[0 0],'linestyle','--','color','k')
line([0 0],[-0.6 0.8],'linestyle','--','color','k')
legend ([plots(1) plots(2) plots(3) plots(4)],titles)
xlabel('Time (ms)')
ylabel('Spearman correlation to model (\rho)')
box off
xlabel('Time (ms)')
ylim([-0.5 1])
ccc
%% fMRI-MEG fusion mine
close all;
clc;
for time=1:length(Times)
    for region=regions_order_to_read
        RDM_MEG=reshape(MEG_RDM(:,:,time),[NoC*NoC 1]);
        RDM_MEG=RDM_MEG(~isnan(RDM_MEG));
        RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
        RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));

        RDM_model_Stim_side=reshape(S_model_periph_cor,[NoC*NoC 1]);
        RDM_model_Stim_side=RDM_model_Stim_side(~isnan(RDM_model_Stim_side));
        
        RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
        RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
        
        RDM_model_Stim_periph=reshape(S_model_decod_cor,[NoC*NoC 1]);
        RDM_model_Stim_periph=RDM_model_Stim_periph(~isnan(RDM_model_Stim_periph));
        
        RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
        RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));
               
        
        Correlations_all_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim RDM_model_Stim_side RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');
        Correlations_all_expt_stim_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim_side RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');
        Correlations_all_expt_stim_side_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');
        Correlations_all_expt_stim_periph_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim_side RDM_model_Stim RDM_model_Rule],'type','Spearman');
        Correlations_all_expt_rule_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim RDM_model_Stim_side RDM_model_Stim_periph],'type','Spearman');
        
        R_sqrd(1,region,time)=Correlations_all_expt_stim_side_out(time,region)-Correlations_all_out(time,region);
        R_sqrd(2,region,time)=Correlations_all_expt_stim_out(time,region)-Correlations_all_out(time,region);
        R_sqrd(3,region,time)=Correlations_all_expt_stim_periph_out(time,region)-Correlations_all_out(time,region);
        R_sqrd(4,region,time)=Correlations_all_expt_rule_out(time,region)-Correlations_all_out(time,region);
        
%         Correlations_all_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim_side RDM_model_Stim_periph RDM_model_Rule],'type','Spearman');
%         Correlations_all_expt_stim_periph_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim_side RDM_model_Rule],'type','Spearman');
%         Correlations_all_expt_rule_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim_periph RDM_model_Stim_side],'type','Spearman');
%         R_sqrd(1,region,time)=Correlations_all_expt_stim_periph_out(time,region)-Correlations_all_out(time,region);
%         R_sqrd(2,region,time)=Correlations_all_expt_rule_out(time,region)-Correlations_all_out(time,region);
    end
end

%% Plotting Correlations my way
close all;
figure;
colors={[0 0 0],[1 0 0],[0.7 0 0],[0 0 1],[0 0 0.7],[1 0 1],[0.7 0 0.7],...
    [1 1 0],[0.7 0.7 0],[0 1 1],[0 0.7 0.7]};
titles={'Stimulus side','Stimulus','Stimulus periph','Rule'};

R_sqrd_summerized=[R_sqrd(:,1,:) nanmean(R_sqrd(:,2:3,:),2) nanmean(R_sqrd(:,4:5,:),2),...
    nanmean(R_sqrd(:,6:7,:),2) nanmean(R_sqrd(:,8:9,:),2) nanmean(R_sqrd(:,10:11,:),2)];

if strcmp(Stim_Resp,'Stim')
    times=[-500:25:3995];
else
    times=[-4000:25:995];
end

c=0;
for region=1:length(Regions_summereized)
    c=c+1;
    plot(times,smooth(squeeze(R_sqrd(4,region,:)),10)','LineWidth',3,'color',colors{c});
    hold on;
end
ylim([-1 1])
line([min(times) max(times)],[0 0],'linestyle','--','color','k')
line([0 0],[-1 1],'linestyle','--','color','k')
xlabel('Time (ms)')
ylabel('Correlation (\rho)')
box off
set(gca,'FontSize',18)
legend(Regions_summereized)


%% fMRI-MEG fusion Martin Hobart
addpath('Z:/projects/Hamid/Projects/MD/Analyses\elife-32816-code1-v2\helper_functions');
disp('Running commonality analysis...')
for time=1:length(Times)
    for region=regions_order_to_read
        RDM_MEG=reshape(MEG_RDM(:,:,time),[NoC*NoC 1]);
        RDM_MEG=RDM_MEG(~isnan(RDM_MEG));
        RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
        RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
        
        RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
        RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));        
        
        RDM_model_Stim_side=reshape(S_model_periph_cor,[NoC*NoC 1]);
        RDM_model_Stim_side=RDM_model_Stim_side(~isnan(RDM_model_Stim_side));
        
        RDM_model_Stim_periph=reshape(S_model_decod_cor,[NoC*NoC 1]);
        RDM_model_Stim_periph=RDM_model_Stim_periph(~isnan(RDM_model_Stim_periph));
        
        
        RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
        RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));
        
        
        xMRI = RDM_fMRI;
        yMEG = RDM_MEG;
        
        
        xstim = RDM_model_Stim;
        xstim_s = RDM_model_Stim_side;
        xstim_p = RDM_model_Stim_periph;
        xrule = RDM_model_Rule;

        
       % % and now again with another calculation
        rMEG_MRIall = correlate([yMEG xMRI xstim_s xstim xstim_p xrule],'type','spearman','method','semipartialcorr');
        rMEG_MRIrule = correlate([yMEG xMRI xstim_s xstim xstim_p],'type','spearman','method','semipartialcorr');
        rMEG_MRIstim_s = correlate([yMEG xMRI xstim xstim_p xrule],'type','spearman','method','semipartialcorr');
        rMEG_MRIstim = correlate([yMEG xMRI xstim_s xstim_p xrule],'type','spearman','method','semipartialcorr');
        rMEG_MRIstim_p = correlate([yMEG xMRI xstim_s xstim xrule],'type','spearman','method','semipartialcorr');
        
        CMEGMRI(1,region,time) = rMEG_MRIstim_s(2,1).^2-rMEG_MRIall(2,1).^2;
        CMEGMRI(2,region,time) = rMEG_MRIstim(2,1).^2-rMEG_MRIall(2,1).^2;
        CMEGMRI(3,region,time) = rMEG_MRIstim_p(2,1).^2-rMEG_MRIall(2,1).^2;
        CMEGMRI(4,region,time) = rMEG_MRIrule(2,1).^2-rMEG_MRIall(2,1).^2;
        
%         rMEG_MRIall = correlate([yMEG xMRI xstim_s xstim_p xrule],'type','spearman','method','semipartialcorr');
%         rMEG_MRIrule = correlate([yMEG xMRI xstim_p xstim_s],'type','spearman','method','semipartialcorr');
%         rMEG_MRIstim_s = correlate([yMEG xMRI xrule xstim_s],'type','spearman','method','semipartialcorr');
%         CMEGMRI(1,region,time) = rMEG_MRIstim_s(2,1).^2-rMEG_MRIall(2,1).^2;
%         CMEGMRI(2,region,time) = rMEG_MRIrule(2,1).^2-rMEG_MRIall(2,1).^2;       
        
    end
end
disp('done.')
%% Plotting commonality
figure;
colors={[0 0 0],[1 0 0],[0.7 0 0],[0 0 1],[0 0 0.7],[1 0 1],[0.7 0 0.7],...
    [1 1 0],[0.7 0.7 0],[0 1 1],[0 0.7 0.7]};

CMEGMRI_summerized=[CMEGMRI(:,1,:) nanmean(CMEGMRI(:,2:3,:),2) nanmean(CMEGMRI(:,4:5,:),2),...
    nanmean(CMEGMRI(:,6:7,:),2) nanmean(CMEGMRI(:,8:9,:),2) nanmean(CMEGMRI(:,10:11,:),2)];

if strcmp(Stim_Resp,'Stim')
    times=[-500:25:3995];
else    
    times=[-4000:25:995];
end
c=0;
for region=1:length(Regions_summereized)
    c=c+1;
    plot(times,smooth(squeeze(CMEGMRI_summerized(1,region,:)),10),'LineWidth',3,'color',colors{c});
    hold on;
end
% ylim([-0.1 0.1])
set(gca,'fontsize',18)
line([min(times) max(times)],[0 0],'linestyle','--','color','k')
line([0 0],[-0.1 0.1],'linestyle','--','color','k')
xlabel('Time (ms)')
ylabel('Commonality coefficient')
box off
set(gca,'FontSize',18)
legend(Regions_summereized)
