clc;
% close all;
clear all;
Main_analysis_directory='D:\Hamid\Postdoc\MD\Analyses';
barwidths=0.3;

subjs=3; %1=all; 2=rule error; 3= stim error; 4=common; 5=one-handed 6= two-handed
excluding_two_handed=0; %1=yes; 0=no

% subjects=[6:26 28:31]; %subjects who did two-handed
subjects=[23 12 19 20 26];


subjects2=[1:5]; %subjects who did two-handed

summerize=1;

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
Regions={'ACC SMA','Left IPS','Right IPS','Left DLPFC','Right DLPFC',...
    'Left VLPFC','Right VLPFC','Left LOC','Right LOC','Left Visual','Right Visual'};

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

    
    Decodings2=nan*ones(31,11);
    for subject=subjects2
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        for region=[1:11]
            load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
            if ~isempty(Decoding_results{region,iteration})
                Decodings2(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    
    data1=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    xtickpositions=[1:size(data1,2)]-barwidths;
    
    data2=Decodings2(:,[1 3 6 2 5 4 7 8 9 10 11]);
    xtickpositions2=[1:size(data2,2)]-barwidths+0.3;

    if summerize==1
        data1=[data1(:,1) nanmean(data1(:,2:3),2) nanmean(data1(:,4:5),2),...
            nanmean(data1(:,6:7),2) nanmean(data1(:,8:9),2) nanmean(data1(:,10:11),2)];
        data2=[data2(:,1) nanmean(data2(:,2:3),2) nanmean(data2(:,4:5),2),...
            nanmean(data2(:,6:7),2) nanmean(data2(:,8:9),2) nanmean(data2(:,10:11),2)];
        xtickpositions=[1:size(data1,2)]-barwidths;
        xtickpositions2=[1:size(data2,2)]-barwidths+0.3;
        Regions={'ACC/pre-SMA','IPS','IFS','AI/FO','LOC','Visual'};
    end
    
    
    subplot(1,2,info);
    bars(1)=bar(xtickpositions,nanmean(data1),'facecolor','g','barwidth',barwidths);
    hold on
    errorbar(xtickpositions,nanmean(data1),nanstd(data1)./sqrt(sum(~isnan(data1(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
   
    bars(2)=bar(xtickpositions2,nanmean(data2),'facecolor','r','barwidth',barwidths);
    hold on
    errorbar(xtickpositions2,nanmean(data2),nanstd(data2)./sqrt(sum(~isnan(data2(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
    
    for i=1:size(data1,2)
         [t,p]=ttest2(data1(~isnan(data1(:,i)),i),data2(~isnan(data2(:,i)),i));
         if t==1
             plot(xtickpositions2(i)-0.15,11,'*b')
         end            
    end
    
    for i=1:size(data1,2)
        data1_temp=data1(~isnan(data1(:,i)),i);
        combinations=nchoosek([1:length(subjects)],5);
        for comb=1:size(combinations,1)
            data1_cmb_mean(comb)=mean(data1_temp(combinations(comb,:)));
        end
        if (nanmean(data2(~isnan(data2(:,i)),i))>0 && sum(data1_cmb_mean<nanmean(data2(~isnan(data2(:,i)),i)))>0.95*length(combinations)) || (nanmean(data2(~isnan(data2(:,i)),i))<0 && sum(data1_cmb_mean>nanmean(data2(~isnan(data2(:,i)),i)))>0.95*length(combinations))
            plot(xtickpositions2(i)-0.15,10,'*m')
        end
    end
    
    xticks([xtickpositions2-0.15])
    xticklabels(Regions)
    xtickangle(45)
    yticks([-20:10:50])
    yticklabels({'30','40','50','60','70','80','90','100'})
    ylim([-20 50])
    ylabel('Accuracy (%)')
    title(titles{info})
    legend ([bars(1) bars(2)],{'One-handed','Two-handed'},'location','northwest')
    box off; grid on;    
end


