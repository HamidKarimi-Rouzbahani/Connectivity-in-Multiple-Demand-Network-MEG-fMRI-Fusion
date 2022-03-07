clc;
close all;
clear all;
Main_analysis_directory='D:/Hamid/Postdoc/MD/Analyses';

Action=4;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs

Trial_type=3; % 1=correct (for RSA=colors*rules*stimuli); 2=same errors (for RSA= the same as previous);
% 3=swapped errors (only for decoding)

% 3=rule coding across cue colors (only for RSA: 1 to 4 cue colors collapsed across all stimuli)
% 4=stimulus coding  (only for RSA: 1 to 4 stimuli collapsed across rules and cues):
% 5=peripheral stimulus coding (only for RSA: collapsed across stimulus side and cue colors from different rules)


Stim_Resp='Stim';
subjects=[1:26 28:31];



addpath('D:\Hamid\Postdoc\MD\MEG')
if Trial_type==3
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_cue_coding_RDM.mat']);
    sampling_time=5;
    
    %     load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_cue_coding_RDM_dec.mat']);
    %     sampling_time=25;

    % blue, green, orange and pink
elseif Trial_type==4
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_stm_coding_RDM.mat']);
    sampling_time=5;
elseif Trial_type==5
    
    load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Stm_coding_Clpsd_cue_colors_RDM.mat']);
    sampling_time=5;
    
    %     load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Stm_coding_Clpsd_cue_colors_RDM_dec.mat']);
    %     sampling_time=25;
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
        if Trial_type==3
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Cue_Color_ROIs.mat']),'RSA_results');
        elseif Trial_type==4
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Stm_ROIs.mat']),'RSA_results');
        elseif Trial_type==5
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['RSA_Stm_clpsdCueCol_ROIs.mat']),'RSA_results');
        end
        %         combinations=nchoosek([1:NoC],2);
        %         for i=1:length(combinations)
        %             if isempty(RSA_results{region, i})
        %                 tmp_RSA(combinations(i,1),combinations(i,2),c)=nan;
        %             else
        %                 tmp_RSA(combinations(i,1),combinations(i,2),c)=RSA_results{region, i}.accuracy.output;
        %             end
        %         end
        tmp_RSA(:,:,c)=1-RSA_results{1, region}.other.output{1,1};
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
load(['D:\Hamid\Postdoc\MD\MEG\RDM_models_',num2str(NoC),'.mat'])
if Trial_type==4
    S_model_cor=nan(NoC);
    S_model_cor(1,4)=1;
    S_model_cor(1,2:3)=0;
    S_model_cor(2,4)=0;
    S_model_cor(3,4)=0;
    S_model_cor(2,3)=1;
end

%% Plotting fMRI and mean MEG RDMs
figure;
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

if Trial_type==3
    subplot_tmp=subplot(4,4,13);
    hold(subplot_tmp,'on');
    image(R_model_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    title('Rule Model')
    
elseif Trial_type>3
    subplot_tmp=subplot(4,4,13);
    hold(subplot_tmp,'on');
    image(S_model_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    title('Stim periph Model')
    title('Stim periph Model')
    
end
%% plotting fMRI correlation to model RDMs
figure;
if Trial_type==3
    titles={'Rule'};
elseif Trial_type>3
    titles={'Stim periph'};
end
c=0;
for region=regions_order_to_read
    c=c+1;
    RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
    RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
    
    if Trial_type>3
        RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
        RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
        cors(1,c)=corr(RDM_fMRI,RDM_model_Stim,'type','Pearson');
    elseif Trial_type==3
        RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
        RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));
        cors(1,c)=corr(RDM_fMRI,RDM_model_Rule,'type','Pearson');
    end
    
end
cors_summerized=[cors(:,1) nanmean(cors(:,2:3),2) nanmean(cors(:,4:5),2),...
    nanmean(cors(:,6:7),2) nanmean(cors(:,8:9),2) nanmean(cors(:,10:11),2)];

subplot(2,1,1);
bar([1:size(cors,2)],cors(1,:))
hold on
xticks([1:length(Regions)])
xticklabels(Regions)
xtickangle(45)
% ylim([-0.8 1])
ylabel('Pearson''s correlation to model(\rho)')
title(titles)
box off

subplot(2,1,2);
bar([1:size(cors_summerized,2)],cors_summerized(1,:))
hold on
xticks([1:length(Regions_summereized)])
xticklabels(Regions_summereized)
xtickangle(45)
% ylim([-0.8 1])
ylabel('Pearson''s correlation to model(\rho)')
title(titles)
box off
%% Plotting MEG correlation to model RDMs
figure;
MEG_RDM=RDM_correct_4_all;
MEG_RDM_err=RDM_error_4_all;

if strcmp(Stim_Resp,'Stim')
    Times=[-500:sampling_time:3995];
else
    Times=[-4000:sampling_time:995];
end

for time=1:length(Times)
    if Trial_type>3
        correl(time)=corr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(S_model_cor,[NoC*NoC 1]),'type','Pearson','rows','complete');
    else
        correl(time)=corr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(R_model_cor,[NoC*NoC 1]),'type','Pearson','rows','complete');
    end
end
if Trial_type==3
    titles={'Rule'};
elseif Trial_type>3
    titles={'Stim periph'};
end

plots=plot(Times,smooth(correl,10),'linewidth',2);

line([min(Times) max(Times)],[0 0],'linestyle','--','color','k')
line([0 0],[-0.6 0.8],'linestyle','--','color','k')
legend (plots,titles)
xlabel('Time (ms)')
ylabel('Spearman correlation to model (\rho)')
box off
xlabel('Time (ms)')

%% fMRI-MEG fusion mine
clc;
c=0;
for region=regions_order_to_read
    c=c+1;
    for time=1:length(Times)
        RDM_MEG=reshape(MEG_RDM(:,:,time),[NoC*NoC 1]);
        RDM_MEG=RDM_MEG(~isnan(RDM_MEG));
        RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
        RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
        
        RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
        RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
        
        RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
        RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));
        Correlations(time,region)=corr(RDM_MEG,RDM_fMRI,'type','Pearson');
        
        data(time,c)=Correlations(time,region);
        
        %         Correlations_all_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim RDM_model_Rule],'type','Pearson');
        %         Correlations_all_expt_rule_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,RDM_model_Stim,'type','Pearson');
        %         Correlations_all_expt_stim_out(time,region)=partialcorr(RDM_MEG,RDM_fMRI,RDM_model_Rule,'type','Pearson');
        %         R_sqrd_stim(time,region)=Correlations_all_expt_stim_out(time,region)-Correlations_all_out(time,region);
        %         R_sqrd_rule(time,region)=Correlations_all_expt_rule_out(time,region)-Correlations_all_out(time,region);
        % %
        % %         R_sqrd_stim2(time,region)=Correlations(time,region)-Correlations_all_expt_stim_out(time,region);
        % %         R_sqrd_rule2(time,region)=Correlations(time,region)-Correlations_all_expt_rule_out(time,region);
    end
end

%% Plotting Correlations my way
figure;
colors={[0 0 0],[1 0 0],[0.7 0 0],[0 0 1],[0 0 0.7],[1 0 1],[0.7 0 0.7],...
    [1 1 0],[0.7 0.7 0],[0 1 1],[0 0.7 0.7]};

if strcmp(Stim_Resp,'Stim')
    times=[-500:sampling_time:3995];
else
    times=[-4000:sampling_time:995];
end
subplot(1,2,1)
c=0;
for region=regions_order_to_read
    c=c+1;
    plot(times,smooth(Correlations(:,region),30),'LineWidth',3,'color',colors{c});
    hold on;
end
% ylim([-1 1])
set(gca,'fontsize',18)
line([min(times) max(times)],[0 0],'linestyle','--','color','k')
line([0 0],[-1 1],'linestyle','--','color','k')
ylabel('Pearson''s correlation (\rho)')
box off
set(gca,'FontSize',18)
legend(Regions)

subplot(1,2,2)
cors_summerized=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
    nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
colors_summerized={[0 0 0],[1 0 0],[0 0 1],[1 0 1],[1 1 0],[0 1 1]};
c=0;
for region=1:size(cors_summerized,2)
    c=c+1;
    plot(times,smooth(cors_summerized(:,region),30),'LineWidth',3,'color',colors_summerized{c});
    hold on;
end
% ylim([-1 1])
set(gca,'fontsize',18)
line([min(times) max(times)],[0 0],'linestyle','--','color','k')
line([0 0],[-1 1],'linestyle','--','color','k')
xlabel('Time (ms)')
ylabel('Pearson''s correlation (\rho)')
box off
legend(Regions_summereized)
ccc

%% fMRI-MEG fusion Martin Hobart
addpath('D:\Hamid\Postdoc\MD\Analyses\elife-32816-code1-v2\helper_functions');
disp('Running commonality analysis...')

for time=1:length(Times)
    c=0;
    for region=regions_order_to_read
        c=c+1;
        RDM_MEG=reshape(MEG_RDM(:,:,time),[NoC*NoC 1]);
        RDM_MEG=RDM_MEG(~isnan(RDM_MEG));
        RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
        RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
        
        RDM_model_Stim=reshape(S_model_cor,[NoC*NoC 1]);
        RDM_model_Stim=RDM_model_Stim(~isnan(RDM_model_Stim));
        
        RDM_model_Rule=reshape(R_model_cor,[NoC*NoC 1]);
        RDM_model_Rule=RDM_model_Rule(~isnan(RDM_model_Rule));
        
        
        xMRI = RDM_fMRI;
        yMEG = RDM_MEG;
        
        
        xstim = RDM_model_Stim;
        xrule = RDM_model_Rule;
        
        if Trial_type==4 || Trial_type==5
            %             rMEG_MRIall = correlate([yMEG xMRI xstim xrule],'type','Pearson','method','semipartialcorr');
            %             rMEG_MRIrule = correlate([yMEG xMRI xrule],'type','Pearson','method','semipartialcorr');
            %             rMEG_MRIstim = correlate([yMEG xMRI xstim],'type','Pearson','method','semipartialcorr');
            %
            %             CMEGMRIstim(time,region) = rMEG_MRIrule(2,1).^2-rMEG_MRIall(2,1).^2;
            %             CMEGMRIrule(time,region) = rMEG_MRIstim(2,1).^2-rMEG_MRIall(2,1).^2;
            
            rMEG_MRIall = correlate([yMEG xMRI xstim],'type','Pearson','method','semipartialcorr');
            rMEG_MRIstim = correlate([yMEG xMRI],'type','Pearson','method','semipartialcorr');
            CMEGMRIstim(time,region) = rMEG_MRIstim(2,1).^2-rMEG_MRIall(2,1).^2;
            data(time,c)=CMEGMRIstim(time,region);
        elseif Trial_type==3
            rMEG_MRIall = correlate([yMEG xMRI xrule],'type','Pearson','method','semipartialcorr');
            rMEG_MRIrule = correlate([yMEG xMRI],'type','Pearson','method','semipartialcorr');
            CMEGMRIrule(time,region) = rMEG_MRIrule(2,1).^2-rMEG_MRIall(2,1).^2;
            data(time,c)=CMEGMRIrule(time,region);
        end
        
    end
end
disp('done.')
%% Plotting commonality
figure;
colors={[0 0 0],[1 0 0],[0.7 0 0],[0 0 1],[0 0 0.7],[1 0 1],[0.7 0 0.7],...
    [1 1 0],[0.7 0.7 0],[0 1 1],[0 0.7 0.7]};

if strcmp(Stim_Resp,'Stim')
    times=[-500:sampling_time:3995];
else
    times=[-4000:sampling_time:995];
end
subplot(1,2,1)
c=0;
for region=regions_order_to_read
    c=c+1;
    if Trial_type>3
        plot(times,smooth(CMEGMRIstim(:,region),30),'LineWidth',3,'color',colors{c});
    elseif Trial_type==3
        plot(times,smooth(CMEGMRIrule(:,region),30),'LineWidth',3,'color',colors{c});
    end
    hold on;
end
set(gca,'fontsize',18)
line([min(times) max(times)],[0 0],'linestyle','--','color','k')
line([0 0],[-0.3 0.8],'linestyle','--','color','k')
ylabel('Commonality coefficient')
box off
set(gca,'FontSize',18)
legend(Regions)

colors_summerized={[0 0 0],[1 0 0],[0 0 1],[1 0 1],[1 1 0],[0 1 1]};

subplot(1,2,2)
if Trial_type>3
    cors_summerized=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
elseif Trial_type==3
    cors_summerized=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
end

c=0;
for region=1:size(cors_summerized,2)
    c=c+1;
    plot(times,smooth(cors_summerized(:,region),30),'LineWidth',3,'color',colors_summerized{c});
    hold on;
end
set(gca,'fontsize',18)
line([min(times) max(times)],[0 0],'linestyle','--','color','k')
line([0 0],[-1 1],'linestyle','--','color','k')
xlabel('Time (ms)')
ylabel('Commonality coefficient')
box off
legend(Regions_summereized)
