clc;
close all;
clear all;
Main_analysis_directory='D:/Hamid/Postdoc/MD/Analyses';

Action=4;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs

Trial_type=11; % 1=correct (for RSA=colors*rules*stimuli); 2=same errors (for RSA= the same as previous);
tt=0;
Stim_Resp='Stim';
iterations=1000;
datas={'StimSide','StimPeriph','Rule','Response'};
for Trial_type=[10 1 3 11]
    tt=tt+1;
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
    
    
    smoothing_rate=100;
    
    
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
        load(['D:\Hamid\Postdoc\MD\MEG\Subjects_Average_',Stim_Resp,'_aligned_Resp_coding_InOut_acrs_opst_col2_RDM_decoding.mat']);
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
    load(['D:\Hamid\Postdoc\MD\MEG\RDM_models_',num2str(NoC),'.mat'])
    if Trial_type==4
        %     % inner vs outer
        %         S_model_cor=nan(NoC);
        %         S_model_cor(1,4)=0;
        %         S_model_cor(1,2:3)=1;
        %         S_model_cor(2,4)=1;
        %         S_model_cor(3,4)=1;
        %         S_model_cor(2,3)=0;
        
        % % side
        S_model_cor=nan(NoC);
        S_model_cor(1,3:4)=1;
        S_model_cor(2,3:4)=1;
        S_model_cor(1,2)=0;
        S_model_cor(3,4)=0;
    end
    if Trial_type==7
        S_model_cor=nan(NoC);
        %     S_model_cor(1,2:3)=0;
        %     S_model_cor(1,4)=1;
        %     S_model_cor(2,3)=1;
        %     S_model_cor(2,4)=0;
        %     S_model_cor(3,4)=0;
        
        S_model_cor(1,2:3)=0;
        S_model_cor(1,4)=1;
        S_model_cor(2,3)=1;
        S_model_cor(2,4)=0;
        S_model_cor(3,4)=0;
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
        
    elseif Trial_type==1 || Trial_type>3
        subplot_tmp=subplot(4,4,13);
        hold(subplot_tmp,'on');
        image(S_model_cor,'Parent',subplot_tmp,'CDataMapping','scaled');
        axis(subplot_tmp,'tight');
        axis(subplot_tmp,'ij');
        set(subplot_tmp,'CLim',scales,'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
        if Trial_type==4
            title('Stim side Model')
        elseif Trial_type==1 || Trial_type==5
            title('Stim periph Model')
        end
        if Trial_type==7
            title('Resp periph Model')
        end
    end
    
    %% plotting fMRI correlation to model RDMs
    figure;
    if Trial_type==3
        titles={'Rule'};
    elseif Trial_type==1 || Trial_type>3
        titles={'Stim periph'};
    end
    c=0;
    for region=regions_order_to_read
        c=c+1;
        RDM_fMRI=reshape(RSA_fMRI(:,:,region),[NoC*NoC 1]);
        RDM_fMRI=RDM_fMRI(~isnan(RDM_fMRI));
        
        if Trial_type==1 || Trial_type>3
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
        Times=[-500:sampling_time:3995]+sampling_time./2;
    else
        Times=[-4000:sampling_time:995]+sampling_time./2;
    end
    
    for time=1:length(Times)
        if Trial_type==1 || Trial_type>3
            correl_MEG(time)=corr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(S_model_cor,[NoC*NoC 1]),'type','Pearson','rows','complete');
        elseif Trial_type==3
            correl_MEG(time)=corr(reshape(MEG_RDM(:,:,time),[NoC*NoC 1]),reshape(R_model_cor,[NoC*NoC 1]),'type','Pearson','rows','complete');
        end
    end
    if Trial_type==3
        titles={'Rule'};
    elseif  Trial_type==1 || Trial_type>3
        titles='Stim periph';
    end
    
    plots=plot(Times,smooth(correl_MEG,smoothing_rate),'linewidth',2);
    
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
            
            if    Trial_type==3
                Correlations_all_out(time,region,1)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Rule],'type','Pearson');
                Correlations_all_expt_target_out(time,region,1)=corr(RDM_MEG,RDM_fMRI,'type','Pearson');
                fusion(time,c,1)=Correlations_all_expt_target_out(time,region,1)-Correlations_all_out(time,region,1);
            elseif Trial_type==1 || Trial_type>3
                Correlations_all_out(time,region,1)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim],'type','Pearson');
                Correlations_all_expt_target_out(time,region,1)=corr(RDM_MEG,RDM_fMRI,'type','Pearson');
                fusion(time,c,1)=Correlations_all_expt_target_out(time,region,1)-Correlations_all_out(time,region,1);
            end
            
            for iteration=2:iterations
                    RDM_model_Rule_shuffled=randsample(RDM_model_Rule,length(RDM_model_Rule));
                    RDM_model_Stim_shuffled=randsample(RDM_model_Stim,length(RDM_model_Stim));
                if    Trial_type==3
                    Correlations_all_out(time,region,iteration)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Rule_shuffled],'type','Pearson');
                    Correlations_all_expt_target_out(time,region,iteration)=corr(RDM_MEG,RDM_fMRI,'type','Pearson');
                    fusion(time,c,iteration)=Correlations_all_expt_target_out(time,region,iteration)-Correlations_all_out(time,region,iteration);
                elseif Trial_type==1 || Trial_type>3
                    Correlations_all_out(time,region,iteration)=partialcorr(RDM_MEG,RDM_fMRI,[RDM_model_Stim_shuffled],'type','Pearson');
                    Correlations_all_expt_target_out(time,region,iteration)=corr(RDM_MEG,RDM_fMRI,'type','Pearson');
                    fusion(time,c,iteration)=Correlations_all_expt_target_out(time,region,iteration)-Correlations_all_out(time,region,iteration);
                end
            end
            [Trial_type region time]
        end
    end
    save(['Fusion_',datas{tt},'_',Stim_Resp,'_aligned.mat'],'Correlations_all_out','Correlations_all_expt_target_out','fusion')
end
ccc
%% Plotting Correlations my way
clc;
close all;
clear all;
Regions_summereized={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};
datas={'StimSide','StimPeriph','Rule','Response'};
iteration=3;
figure;
Stim_Resp='Resp';
sampling_time=5;
smoothing_rate=100;
if strcmp(Stim_Resp,'Stim')
    Times=[-500:sampling_time:3995]+(sampling_time*5)./2+(smoothing_rate*5)./2;
else
    Times=[-4000:sampling_time:995]+(sampling_time*5)./2-(smoothing_rate*5)./2;
end
colors_summerized={[1 0 0],[0 1 0],[0 0 1],[0 0 0]};
for tt=1:4
    load(['Fusion_',datas{tt},'_',Stim_Resp,'_aligned.mat']);
    cors_summerizedNew(:,:,tt)=[fusion(:,1,iteration) nanmean(fusion(:,2:3,iteration),2) nanmean(fusion(:,4:5,iteration),2),...
        nanmean(fusion(:,6:7,iteration),2) nanmean(fusion(:,8:9,iteration),2) nanmean(fusion(:,10:11,iteration),2)];
end

legends={'Stim Side','Stim Periph','Rule','Resp'};
c=0;
regions_order_to_plot=[4 3 1 6 5 2];
for region=regions_order_to_plot
    c=c+1;
    subplot(2,3,c)
    for dtt=1:4
        plots(dtt)=plot(Times,smooth(cors_summerizedNew(:,region,dtt),smoothing_rate),'LineWidth',3,'color',colors_summerized{dtt});
        hold on;
        ylim([-1 1.3])
        set(gca,'fontsize',18)
        line([min(Times) max(Times)],[0 0],'linestyle','--','color','k')
        line([0 0],[-1 1.3],'linestyle','--','color','k')
        if strcmp(Stim_Resp,'Stim')
            xlim([-200 3500])
        else
            xlim([-3500 500])
        end
        box off
    end
    if c==1
        legend([plots(1) plots(2) plots(3) plots(4)],legends,'location','northeast','fontsize',10)
    end
    if c==1 || c==4
        ylabel('Pearson''s correlation (r)')
    end
    if c>3
        xlabel('Time (ms)')
    end
    title(Regions_summereized{region})
end

