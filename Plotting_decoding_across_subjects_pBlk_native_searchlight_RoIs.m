clc;
% close all;
clear all;
Main_analysis_directory='D:\Hamid\Postdoc\MD\Analyses';
barwidths=0.3;

Action=1;
Trial_type=1;

subjects=[1:26 28:31];
% subjects=[1:5];


if Action==1
    Beta_folder_name='Correct_trials_decoding_design_pBlk_native';
elseif Action==2
    Beta_folder_name='Correct_trials_RSA_design_16_pBlk_native';
    categories=16;
elseif Action==3
    Beta_folder_name='Correct_trials_RSA_design_8_pBlk_native';
    categories=8;
elseif Action==4
    Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
    categories=4;
end
if Trial_type==2
    Beta_folder_name=[Beta_folder_name,'_AllErr'];
elseif Trial_type==3
    Beta_folder_name=[Beta_folder_name,'_OpsErr'];
end

Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];

Regions={'Rule region1','Rule region2','Stimulus region1','Stimulus region2'};


Action=1;
Trial_type=1;

aspects={'S','R'};
titles={'Stimulus decoding','Rule decoding'};
iteration=1;
for info=1:2
    Decodings=nan*ones(31,4);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:4]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_my_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings;
    xtickpositions=[1:size(data,2)];
    
    if info==2
        decoding_averaged_for_correlation=data;
    end
    
    subplot(1,2,info);
    bars(1)=bar(xtickpositions,nanmean(data),'facecolor','g','barwidth',barwidths);
    hold on
    errorbar(xtickpositions,nanmean(data),nanstd(data)./sqrt(sum(~isnan(data(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions)])
    xticklabels(Regions)
    xtickangle(45)
    yticks([-20:10:50])
    yticklabels({'30','40','50','60','70','80','90','100'})
    ylim([-20 50])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off; grid on;
end

colors={'r','b'};
Action=1;
Trial_type=2;
if Action==1
    Beta_folder_name='Correct_trials_decoding_design_pBlk_native';
elseif Action==2
    Beta_folder_name='Correct_trials_RSA_design_16_pBlk_native';
    categories=16;
elseif Action==3
    Beta_folder_name='Correct_trials_RSA_design_8_pBlk_native';
    categories=8;
elseif Action==4
    Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
    categories=4;
end
if Trial_type==2
    Beta_folder_name=[Beta_folder_name,'_AllErr'];
elseif Trial_type==3
    Beta_folder_name=[Beta_folder_name,'_OpsErr'];
end
aspects={'S','R'};
titles={'Stimulus decoding','Rule decoding'};
iteration=1;
for info=1:2
    Decodings=nan*ones(31,4);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:4]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_my_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings;
    xtickpositions=[1:size(data,2)]+0.3;
    
    subplot(1,2,info);
    bars(2)=bar(xtickpositions,nanmean(data),'facecolor',colors{info},'barwidth',barwidths);
    hold on
    errorbar(xtickpositions,nanmean(data),nanstd(data)./sqrt(sum(~isnan(data(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions)])
    xticklabels(Regions)
    xtickangle(45)
    yticks([-20:10:50])
    yticklabels({'30','40','50','60','70','80','90','100'})
    ylim([-20 50])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off; grid on;
end


Action=1;
Trial_type=3;
if Action==1
    Beta_folder_name='Correct_trials_decoding_design_pBlk_native';
elseif Action==2
    Beta_folder_name='Correct_trials_RSA_design_16_pBlk_native';
    categories=16;
elseif Action==3
    Beta_folder_name='Correct_trials_RSA_design_8_pBlk_native';
    categories=8;
elseif Action==4
    Beta_folder_name='Correct_trials_RSA_design_4_pBlk_native';
    categories=4;
end
if Trial_type==2
    Beta_folder_name=[Beta_folder_name,'_AllErr'];
elseif Trial_type==3
    Beta_folder_name=[Beta_folder_name,'_OpsErr'];
end
aspects={'S','R'};
titles={'Stimulus decoding','Rule decoding'};
iteration=1;
for info=1:2
    Decodings=nan*ones(31,4);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:4]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_my_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings;
    xtickpositions=[1:size(data,2)]-0.3;
    
    subplot(1,2,info);
    bars(3)=bar(xtickpositions,nanmean(data),'facecolor',colors{3-info},'barwidth',barwidths);
    hold on
    errorbar(xtickpositions,nanmean(data),nanstd(data)./sqrt(sum(~isnan(data(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions)])
    xticklabels(Regions)
    xtickangle(45)
    yticks([-20:10:50])
    yticklabels({'30','40','50','60','70','80','90','100'})
    ylim([-20 50])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off; grid on;
end

% set('gca','fontsize',18)
legend ([bars(1),bars(2) bars(3)],{'Correct','Rule errors','Stimulus errors'})

ccc
%% correlation to behaviour
clc;
Regions={'Rule region1','Rule region2','Stimulus region1','Stimulus region2'};

load('behavioural_results.mat')
% Thing=Stims_no/0.4;
Thing=Stims_rt;
trial_outcome=[4];
subjects=[1:26 28:31];
for region=1:4
    [r,p]=corr(decoding_averaged_for_correlation(subjects,region),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,subjects)),3),2)),'rows','complete','type','Spearman');
    subplot_tmp=subplot(2,2,region);        
    f=polyfit(decoding_averaged_for_correlation(subjects,region),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,subjects)),3),2)),1);
    scatter(50+decoding_averaged_for_correlation(subjects,region),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,subjects)),3),2)),80,'marker','o','MarkerFaceColor','k');
    set(gca,'fontsize',11);    
    yest=polyval(f,decoding_averaged_for_correlation(subjects,region));
    hold on;
    plot(50+decoding_averaged_for_correlation(subjects,region),yest,'k--','LineWidth',3);    
    hold off;
    xlabel(['Decoding accuracy (%)'])


    title(Regions{region})
    set(subplot_tmp,'FontSize',14,'LineWidth',2)
    xlim([40 100])

% %     ylabel('Behavioral accuracy (%)')
% %     ylim([60 100])
% %     text(87,97,['r = ',num2str(r,'%.2f')],'FontSize',14)
% %     text(87,94,['p = ',num2str(p,'%.2f')],'FontSize',14)
% %     
    ylabel('Reaction time (s)')
    ylim([1 3])
    text(85,2.8,['r = ',num2str(r,'%.2f')],'FontSize',14)
    text(85,2.55,['p = ',num2str(p,'%.2f')],'FontSize',14)
end
% figure;
% scatter(data_averaged_for_correlation(:,region),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,:)),3),2)))

