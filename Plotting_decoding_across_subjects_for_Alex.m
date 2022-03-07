clc;
% close all;
clear all;
Main_analysis_directory='D:/Postdocs/Macqaurie/MD_project/fMRI/Analyses';
barwidths=0.3;

Action=1;
Trial_type=1;
subjects=[1:24 26 28:31];

    if Action==1
        Beta_folder_name='Correct_trials_decoding_design';
    elseif Action==2
        Beta_folder_name='Correct_trials_RSA_design_16';
        categories=16;
    elseif Action==3
        Beta_folder_name='Correct_trials_RSA_design_8';
        categories=8;
    elseif Action==4
        Beta_folder_name='Correct_trials_RSA_design_4';
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
%7,8,16, 25 undone
% subjects=[1:6 9:15 17:22 28:31];

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
    yticks([0:10:50])
    yticklabels({'50','60','70','80','90','100'})
    ylim([-10 35])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off
    
    subplot(2,2,info+2);
    bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor','g','barwidth',barwidths)
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
    xtickangle(45)
    yticks([0:10:50])
    yticklabels({'50','60','70','80','90','100'})
    ylim([-10 35])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off
end

Action=1;
Trial_type=2;
subjects=[1:24 26 28:31];
colors={'r','b'};

    if Action==1
        Beta_folder_name='Correct_trials_decoding_design';
    elseif Action==2
        Beta_folder_name='Correct_trials_RSA_design_16';
        categories=16;
    elseif Action==3
        Beta_folder_name='Correct_trials_RSA_design_8';
        categories=8;
    elseif Action==4
        Beta_folder_name='Correct_trials_RSA_design_4';
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
%7,8,16, 25 undone
% subjects=[1:6 9:15 17:22 28:31];

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
    yticks([0:10:50])
    yticklabels({'50','60','70','80','90','100'})
    ylim([-10 35])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off
    
    subplot(2,2,info+2);
    bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor',colors{info},'barwidth',barwidths)
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
    xtickangle(45)
    yticks([0:10:50])
    yticklabels({'50','60','70','80','90','100'})
    ylim([-10 35])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off
end



Action=1;
Trial_type=3;
subjects=[1:24 26 28];
colors={'r','b'};
    if Action==1
        Beta_folder_name='Correct_trials_decoding_design';
    elseif Action==2
        Beta_folder_name='Correct_trials_RSA_design_16';
        categories=16;
    elseif Action==3
        Beta_folder_name='Correct_trials_RSA_design_8';
        categories=8;
    elseif Action==4
        Beta_folder_name='Correct_trials_RSA_design_4';
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
%7,8,16, 25 undone
% subjects=[1:6 9:15 17:22 28:31];

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
    yticks([0:10:50])
    yticklabels({'50','60','70','80','90','100'})
    ylim([-10 35])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off
    
    subplot(2,2,info+2);
    bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor',colors{3-info},'barwidth',barwidths)
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    xticks([1:length(Regions_summereized)])
    xticklabels(Regions_summereized)
    xtickangle(45)
    yticks([0:10:50])
    yticklabels({'50','60','70','80','90','100'})
    ylim([-10 35])
    ylabel('Accuracy (%)')
    title(titles{info})
    box off
end
    legend([bars(1) bars(2) bars(3)],{'Correct','Rule errors','Stimulus errors'})
ccc
%% Cross-condition Significance Matrix
figure;
for info=1:3
    Decodings=nan*ones(31,11);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:11]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
            Decodings(subject,region)=Decoding_results{1, region}.accuracy_minus_chance.output;
        end
    end
    data=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    Significanc_mat=nan*ones(size(data,2));
    Bayes_mat=nan*ones(size(data,2));
    
    for region1=1:size(data,2)
        for region2=1:size(data,2)
            Significanc_mat(region1,region2)=bf.ttest(squeeze(data(:,region1)),squeeze(data(:,region2)));
            if region1==region2
                Significanc_mat(region1,region2)=1;
            end
        end
    end
    for region1=1:size(data,2)
        for region2=1:size(data,2)
            if Significanc_mat(region1,region2)>10
                Bayes_mat(region1,region2)=6;
            elseif Significanc_mat(region1,region2)>3 && Significanc_mat(region1,region2)<=10
                Bayes_mat(region1,region2)=5;
            elseif Significanc_mat(region1,region2)>1 && Significanc_mat(region1,region2)<=3
                Bayes_mat(region1,region2)=4;
            elseif Significanc_mat(region1,region2)<1 && Significanc_mat(region1,region2)>=1/3
                Bayes_mat(region1,region2)=3;
            elseif Significanc_mat(region1,region2)<1/3 && Significanc_mat(region1,region2)>=1/10
                Bayes_mat(region1,region2)=2;
            elseif Significanc_mat(region1,region2)<1/10
                Bayes_mat(region1,region2)=1;
            end
        end
    end
    subplot_tmp=subplot(2,3,info);
    hold(subplot_tmp,'on');
    image(Bayes_mat(:,:),'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',[1 6],'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    
    xticks(1:size(data,2))
    yticks(1:size(data,2))
    xticklabels(Regions)
    yticklabels(Regions)
    ytickangle(45)
    xtickangle(45)
    colormap(parula(6));
    c = colorbar ('YTickLabel',{'BF<0.1','0.1<BF<0.3','0.3<BF<1','1<BF<3','3<BF<10','BF>10'}) ; %Create Colorbar
    c.Ticks = [1+0.4 2.2 3.1 3.9 4.7 5.6]; %Create ticks
    c.Label.String = 'Bayes Factors';
       title(titles{info})
 
    
     data_averaged=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
  
    Significanc_mat=nan*ones(size(data_averaged,2));
    Bayes_mat=nan*ones(size(data_averaged,2));
    
    for region1=1:size(data_averaged,2)
        for region2=1:size(data_averaged,2)
            Significanc_mat(region1,region2)=bf.ttest(squeeze(data_averaged(:,region1)),squeeze(data_averaged(:,region2)));
            if region1==region2
                Significanc_mat(region1,region2)=1;
            end
        end
    end
    for region1=1:size(data_averaged,2)
        for region2=1:size(data_averaged,2)
            if Significanc_mat(region1,region2)>10
                Bayes_mat(region1,region2)=6;
            elseif Significanc_mat(region1,region2)>3 && Significanc_mat(region1,region2)<=10
                Bayes_mat(region1,region2)=5;
            elseif Significanc_mat(region1,region2)>1 && Significanc_mat(region1,region2)<=3
                Bayes_mat(region1,region2)=4;
            elseif Significanc_mat(region1,region2)<1 && Significanc_mat(region1,region2)>=1/3
                Bayes_mat(region1,region2)=3;
            elseif Significanc_mat(region1,region2)<1/3 && Significanc_mat(region1,region2)>=1/10
                Bayes_mat(region1,region2)=2;
            elseif Significanc_mat(region1,region2)<1/10
                Bayes_mat(region1,region2)=1;
            end
        end
    end
    subplot_tmp=subplot(2,3,info+3);
    hold(subplot_tmp,'on');
    image(Bayes_mat(:,:),'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',[1 6],'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    
    xticks(1:size(data_averaged,2))
    yticks(1:size(data_averaged,2))
    xticklabels(Regions_summereized)
    yticklabels(Regions_summereized)
    ytickangle(45)
    xtickangle(45)
    colormap(parula(6));
    c = colorbar ('YTickLabel',{'BF<0.1','0.1<BF<0.3','0.3<BF<1','1<BF<3','3<BF<10','BF>10'}) ; %Create Colorbar
    c.Ticks = [1+0.4 2.2 3.1 3.9 4.7 5.6]; %Create ticks
    c.Label.String = 'Bayes Factors';   
       title(titles{info}) 
end