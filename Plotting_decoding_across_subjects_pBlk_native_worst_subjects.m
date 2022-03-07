clc;
% close all;
clear all;
Main_analysis_directory='D:\Hamid\Postdoc\MD\Analyses';
barwidths=0.3;

subjs=1; %1=all; 2=rule error; 3= stim error; 4=common; 5=one-handed 6= two-handed 7 = most errors
excluding_two_handed=0; %1=yes; 0=no
   
if subjs==1
    subjects=[1:26 28:31];
elseif subjs==2
    if excluding_two_handed==0
        subjects=[16 17 24 12 10 20 21 2 29 14 23 7 5 26 8]; %rule error subjects
    else
        subjects=[16 17 24 12 10 20 21 29 14 23 7 26 8 25 22]; %rule error subjects
    end
elseif subjs==3
    if excluding_two_handed==0
        subjects=[4 23 12 19 1 20 26 6 10 8 11 15 21 5 3]; %stim error subjects
    else
        subjects=[23 12 19 20 26 6 10 8 11 15 21 22 18 31 17]; %stim error subjects without two-handed
    end
elseif subjs==4
    subjects=[23 20 26 10 8 21 5]; %stim and rule error subjects
elseif subjs==5
    subjects=[6:26 28:31]; %subjects who did one-handed
elseif subjs==6
    subjects=[1:5]; %subjects who did two-handed
elseif subjs==7
    subjects=[14 24 8 5 21 19 1 10 17 20 23 26 16 4 12]; %subjects with most errors
end



Action=1;
Trial_type=1;

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
% Regions={'ACC SMA','Left DLPFC','Left IPS','Left VLPFC',...
%     'Right DLPFC','Right IPS','Right VLPFC'...
%     'Left LOC','Right LOC','Left Visual','Right Visual'};

Regions={'ACC SMA','Left IPS','Right IPS','Left DLPFC','Right DLPFC',...
    'Left VLPFC','Right VLPFC','Left LOC','Right LOC','Left Visual','Right Visual'};

Regions_summereized={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};


aspects={'S','R','P'};
titles={'Stimulus decoding','Rule decoding','Response decoding'};
iteration=1;
for info=1:2
    Decodings=nan*ones(31,11);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:11]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    
    data_averaged=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
    xtickpositions=[1:size(data,2)]-barwidths;
    xtickpositions_smrzd=[1:size(data_averaged,2)]-barwidths;
    
    subplot(2,2,info);
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
    
    subplot(2,2,info+2);
    bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor','g','barwidth',barwidths)
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
    xtickangle(45)
    yticks([-20:10:50])
    yticklabels({'30','40','50','60','70','80','90','100'})
    ylim([-20 50])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off; grid on;
end



Action=1;
Trial_type=2;

if subjs==1 
    subjects=[1:26 28:31];
elseif subjs==2
    if excluding_two_handed==0
        subjects=[16 17 24 12 10 20 21 2 29 14 23 7 5 26 8]; %rule error subjects
    else
        subjects=[16 17 24 12 10 20 21 29 14 23 7 26 8 25 22]; %rule error subjects
    end
elseif subjs==3
    if excluding_two_handed==0
        subjects=[4 23 12 19 1 20 26 6 10 8 11 15 21 5 3]; %stim error subjects
    else
        subjects=[23 12 19 20 26 6 10 8 11 15 21 22 18 31 17]; %stim error subjects without two-handed
    end
elseif subjs==4
    subjects=[23 20 26 10 8 21 5]; %stim and rule error subjects
elseif subjs==5
    subjects=[6:26 28:31]; %subjects who did one-handed
elseif subjs==6
    subjects=[1:5]; %subjects who did two-handed
elseif subjs==7
    subjects=[14 24 8 5 21 19 1 10 17 20 23 26 16 4 12]; %subjects with most errors
end

colors={'r','b'};

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
% Regions={'ACC SMA','Left DLPFC','Left IPS','Left VLPFC',...
%     'Right DLPFC','Right IPS','Right VLPFC'...
%     'Left LOC','Right LOC','Left Visual','Right Visual'};

Regions={'ACC SMA','Left IPS','Right IPS','Left DLPFC','Right DLPFC',...
    'Left VLPFC','Right VLPFC','Left LOC','Right LOC','Left Visual','Right Visual'};

Regions_summereized={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};


aspects={'S','R','P'};
titles={'Stimulus decoding','Rule decoding','Response decoding'};
iteration=1;
for info=1:2
    Decodings=nan*ones(31,11);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:11]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    data_averaged=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
    

    xtickpositions=[1:size(data,2)];
    xtickpositions_smrzd=[1:size(data_averaged,2)];
    
    subplot(2,2,info);
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
    
    subplot(2,2,info+2);
    bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor',colors{info},'barwidth',barwidths)
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
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


if subjs==1 
    subjects=[1:26 28:31];
elseif subjs==2
    if excluding_two_handed==0
        subjects=[16 17 24 12 10 20 21 2 29 14 23 7 5 26 8]; %rule error subjects
    else
        subjects=[16 17 24 12 10 20 21 29 14 23 7 26 8 25 22]; %rule error subjects
    end
elseif subjs==3
    if excluding_two_handed==0
        subjects=[4 23 12 19 1 20 26 6 10 8 11 15 21 5 3]; %stim error subjects
    else
        subjects=[23 12 19 20 26 6 10 8 11 15 21 22 18 31 17]; %stim error subjects without two-handed
    end
elseif subjs==4
    subjects=[23 20 26 10 8 21 5]; %stim and rule error subjects
elseif subjs==5
    subjects=[6:26 28:31]; %subjects who did one-handed
elseif subjs==6
    subjects=[1:5]; %subjects who did two-handed
elseif subjs==7
    subjects=[14 24 8 5 21 19 1 10 17 20 23 26 16 4 12]; %subjects with most errors
end

colors={'r','b'};
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
% Regions={'ACC SMA','Left DLPFC','Left IPS','Left VLPFC',...
%     'Right DLPFC','Right IPS','Right VLPFC'...
%     'Left LOC','Right LOC','Left Visual','Right Visual'};

Regions={'ACC SMA','Left IPS','Right IPS','Left DLPFC','Right DLPFC',...
    'Left VLPFC','Right VLPFC','Left LOC','Right LOC','Left Visual','Right Visual'};

Regions_summereized={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};

aspects={'S','R','P'};
titles={'Stimulus decoding','Rule decoding','Response decoding'};
iteration=1;
for info=1:2
    Decodings=nan*ones(31,11);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:11]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    data_averaged=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
    xtickpositions=[1:size(data,2)]+barwidths;
    xtickpositions_smrzd=[1:size(data_averaged,2)]+barwidths;
    
    if info==2
        data_averaged_for_correlation=data_averaged;
    end
    
    subplot(2,2,info);
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
    
    subplot(2,2,info+2);
    bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor',colors{3-info},'barwidth',barwidths)
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
    xtickangle(45)
    yticks([-20:10:50])
    yticklabels({'30','40','50','60','70','80','90','100'})
    ylim([-20 50])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off; grid on;
end
legend([bars(1) bars(2) bars(3)],{'Correct','Rule errors','Stimulus errors'})
%% correlation to behaviour
clc;
Regions_summereized={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};

load('behavioural_results.mat')
Thing=Stims_no/0.4;
% Thing=Stims_rt;
trial_outcome=[4];
subjects=[6:26 28:31];
for region=1:6
    [r,p]=corr(data_averaged_for_correlation(subjects,region),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,subjects)),3),2)),'rows','complete');
    if p<0.05
       [region p r]
    end
end
% figure;
% scatter(data_averaged_for_correlation(:,region),squeeze(nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,:)),3),2)))

