clc;
% close all;
clear all;
Main_analysis_directory='D:\Hamid\Postdoc\MD\Analyses';
barwidths=0.3;

Action=1;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs

Trial_type=4; % 1=correct (for RSA=colors*rules*stimuli); 2=same errors (for RSA= the same as previous);
% 3=swapped errors (only for decoding) % 4=correct (for decoding=colors*rules*stimuli*stim side)


subjects=[1:26 28:31];

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
elseif Trial_type==4
    Beta_folder_name=[Beta_folder_name,'_IncldSide'];
end

Directory_for_working=[Main_analysis_directory,'/Results_temp_playing/'];
% Regions={'ACC SMA','Left DLPFC','Left IPS','Left VLPFC',...
%     'Right DLPFC','Right IPS','Right VLPFC'...
%     'Left LOC','Right LOC','Left Visual','Right Visual'};

Regions={'ACC SMA','Left IPS','Right IPS','Left DLPFC','Right DLPFC',...
    'Left VLPFC','Right VLPFC','Left LOC','Right LOC','Left Visual','Right Visual'};


Regions_summereized={'Visual','LOC','IPS','AI/FO','IFS','ACC/pre-SMA'};
regions_order_to_plot=[6 5 2 4 3 1];

aspects={'L','S','R','P'};
colors_summerized={[1 0 0],[0 1 0],[0 0 1],[0 0 0]};
iteration=1;

for info=1:4
    Decodings=nan*ones(31,11);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
        for region=[1:11]
            if ~isempty(Decoding_results{region,iteration})
                Decodings(subject,region)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
            end
        end
    end
    data=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    
    data_averaged=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
    data_averaged=data_averaged(:,regions_order_to_plot);
    
    xtickpositions_smrzd=[1:barwidths*10:size(data_averaged,2)*10*barwidths]+info*2*barwidths;
    %     xtickpositions_smrzd=[size(data_averaged,2):-1:1]+info*0.1;
    set(gca,'fontsize',18)
    bars(info)=bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor',colors_summerized{info},'barwidth',barwidths*0.66);
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))),'linewidth',2,'color','k','CapSize',0,'LineStyle','none')
end
xticks([1:barwidths*10:size(data_averaged,2)*10*barwidths]+2.5*2*barwidths)
xticklabels(Regions_summereized)
xtickangle(45)
yticks([-10:10:60])
yticklabels({'','50','60','70','80','90','100',''})
ylim([-5 55])
ylabel('Accuracy (%)')
box off; grid on;
legend([bars(1) bars(2) bars(3) bars(4)],{'Stim side','Stim periph','Rule','Resp'})

%% MEG
clc;
clear all;
close all;
OnsetResponse=1;

for OnsetResponse=1:2
    subplot(1,2,OnsetResponse)
    if OnsetResponse==1
        times=[-500:5:3995];
    elseif OnsetResponse==2
        times=[-4000:5:995];
    end
    colors_summerized={[1 0 0],[0 1 0],[0 0 1],[0 0 0]};
    smoothing_rate=20;
    infos={'Stim_side','Stim_periph','Rule','Resp'};
    set(gca,'fontsize',18)
    for info=1:4
        for Subject=1:19
            if OnsetResponse==1
                load(['p',num2str(Subject),'_Stim_aligned_',infos{info},'_decoding_.mat'],'accuracy_correct');
            elseif OnsetResponse==2
                load(['p',num2str(Subject),'_Resp_aligned_',infos{info},'_decoding.mat'],'accuracy_correct');
            end
            Accuracies(:,Subject)=smooth(accuracy_correct,smoothing_rate);
        end
        plot_line{info}=shadedErrorBar(times,nanmean(Accuracies,2),nanstd(Accuracies')./sqrt(size(Accuracies,2)),{'color',colors_summerized{info},'LineWidth',3},1);
        hold on;
    end
    line([min(times) max(times)],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
    line([0 0],[0.48 0.72],'LineWidth',1.5,'Color','k','LineStyle','--');
    xlim([min(times) max(times)])
    ylim([0.48 0.72])
    
    ylabel('Accuracy (%)')
    box off;
    if OnsetResponse==1
        xlabel('Time Relative to Stimulus Onset (ms)')
        legend([plot_line{1,1}.mainLine plot_line{1,2}.mainLine plot_line{1,3}.mainLine plot_line{1,4}.mainLine],{'Stim side','Stim periph','Rule','Resp'})
    else
        xlabel('Time Relative to Response (ms)')
    end
    box off;
    set(gca,'FontSize',24,'LineWidth',4,'YTick',...
    [0.5 0.6 0.7],'YTickLabel',{'50','60','70'},'XMinorTick','on');
    clearvars Accuracies
end