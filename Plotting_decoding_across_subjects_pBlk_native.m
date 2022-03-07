clc;
% close all;
clear all;
Main_analysis_directory='D:\Hamid\Postdoc\MD\Analyses';
barwidths=0.3;

Action=1;
Trial_type=1;

subjects=[1:26 28:31];
% subjects=[12 16 17 20 24]; %rule error subjects
% subjects=[1 4 12 19 20 23]; %stim error subjects
% subjects=[12 20]; %stim and rule error subjects


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

subjects=[1:26 28:31];
% subjects=[12 16 17 20 24]; %rule error subjects
% subjects=[1 4 12 19 20 23]; %stim error subjects
% subjects=[12 20]; %stim and rule error subjects


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
    bars(2)=bar(xtickpositions,nanmean(data),'facecolor',colors{info},'barwidth',barwidths)
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
subjects=[1:26 28:31];
% subjects=[12 16 17 20 24]; %rule error subjects
% subjects=[1 4 12 19 20 23]; %stim error subjects
% subjects=[12 20]; %stim and rule error subjects

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
