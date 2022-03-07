clc;
% close all;
clear all;
Main_analysis_directory='Z:\projects\Hamid\Projects\MD\Analyses';
barwidths=0.3;
upper_threshold=3;
lower_threshold=upper_threshold./10;

analysis_figures_dir  = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\Figures\New\'];

legend_box_outline     = 'off';
legend_box_loaction    = 'northeast';
legend_fontangel       = 'normal';
legend_fontsize        = 14;
pdf_paper_size=[20 20];
pdf_print_resolution   = '-r300';

Action=1;       % =1 for decoding of stim, rule and resp each two class
% 2:16*16,   3:8*8    4:4*4 RDMs

Trial_type=4; % 1=correct (for RSA=colors*rules*stimuli); 2=same errors (for RSA= the same as previous);
% 3=swapped errors (only for decoding) % 4=correct (for decoding=colors*rules*stimuli*stim side)


subjects=[1:26 28:31];
% subjects=[1:4 17:23];

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


Regions_summereized={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
regions_order_to_plot=[6 5 2 4 3 1];

aspects={'L','S','R','P'};
colors_summerized={[0.5 0.1 0.1],[0.1 0.5 0.1],[0.1 0.1 0.5],[0.3 0.3 0.3]};

data_for_Bayes=nan(4,31,6);
for info=1:4
    Decodings=nan(31,11);
    Decodings_rands=nan(31,11,100);
    for subject=subjects
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        load(fullfile(Directory_for_working,num2str(csub,'S%03d'),'/',Beta_folder_name,'/',['Dec_ROIs_',aspects{info},'.mat']));
        for region=[1:11]
            if ~isempty(Decoding_results{region,1})
                Decodings(subject,region)=Decoding_results{region,1}.accuracy_minus_chance.output;
                iterations=2:101;
                c=0;
                for iteration=iterations
                    c=c+1;
                    Decodings_rands(subject,region,c)=Decoding_results{region,iteration}.accuracy_minus_chance.output;
                end
            end
        end
    end
    data=Decodings(:,[1 3 6 2 5 4 7 8 9 10 11]);
    
    data_averaged=[data(:,1) nanmean(data(:,2:3),2) nanmean(data(:,4:5),2),...
        nanmean(data(:,6:7),2) nanmean(data(:,8:9),2) nanmean(data(:,10:11),2)];
    data_averaged=data_averaged(:,regions_order_to_plot);
    
    data_rand=Decodings_rands(:,[1 3 6 2 5 4 7 8 9 10 11],:);    
    data_averaged_rand=[data_rand(:,1,:) nanmean(data_rand(:,2:3,:),2) nanmean(data_rand(:,4:5,:),2),...
        nanmean(data_rand(:,6:7,:),2) nanmean(data_rand(:,8:9,:),2) nanmean(data_rand(:,10:11,:),2)];
    data_averaged_rand=data_averaged_rand(:,regions_order_to_plot,:);
    
      
    
    xtickpositions_smrzd=[1:barwidths*10:size(data_averaged,2)*10*barwidths]+info*2*barwidths;
    bars(info)=bar(xtickpositions_smrzd,nanmean(data_averaged),'facecolor',colors_summerized{info},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on
    errorbar(xtickpositions_smrzd,nanmean(data_averaged),(nanstd(data_averaged)./sqrt(sum(~isnan(data_averaged(:,1)))))*1.96,'linewidth',2,'color','k','CapSize',0,'LineStyle','none')    
%     plot(xtickpositions_smrzd,significant(info,:)*(-5),'*','color',colors_summerized{info})

    for region=1:6
        significance_bayes_info(region)=bf.ttest2(data_averaged(:,region),nanmean(data_averaged_rand(:,region,:),3));
        num_to_plot=round(significance_bayes_info(region),2,'significant');
        if num_to_plot>upper_threshold || num_to_plot<lower_threshold
            text(xtickpositions_smrzd(region),-14,num2str(num_to_plot,'%10.1e\n'),'color',colors_summerized{info},'Rotation',90,'FontWeight','bold','Fontangle','italic');
        else
            text(xtickpositions_smrzd(region),-14,num2str(num_to_plot,'%10.1e\n'),'color',colors_summerized{info},'Rotation',90,'FontWeight','light','Fontangle','normal');            
        end
    end
    data_for_Bayes(info,:,:)=data_averaged;
end

xtickpositions_smrzd=[1:barwidths*10:size(data_averaged,2)*10*barwidths]+1.4*2*barwidths;
for region=1:6
    significance_bayes_side_periph(region)=bf.ttest2(data_for_Bayes(1,:,region),data_for_Bayes(2,:,region));
    num_to_plot=round(significance_bayes_side_periph(region),2,'significant');
    if num_to_plot>upper_threshold || num_to_plot<lower_threshold
        text(xtickpositions_smrzd(region),52,num2str(num_to_plot,'%10.1e\n'),'color','k','Rotation',90,'FontWeight','bold','Fontangle','italic');
    else
        text(xtickpositions_smrzd(region),52,num2str(num_to_plot,'%10.1e\n'),'color','k','Rotation',90,'FontWeight','light','Fontangle','normal');
    end
end

xticks([1:barwidths*10:size(data_averaged,2)*10*barwidths]+2.5*2*barwidths)
xticklabels(Regions_summereized)
xtickangle(45)
yticks([0:10:50])
yticklabels({'50','60','70','80','90','100'})
ylim([-15 60])
ylabel('Decoding accuracy (%)') 
box off;
plot([0.4 19.6],[0 0],'linewidth',2,'color','k')
% legend([bars(1) bars(2) bars(3) bars(4)],{'Stim side','Stim periph','Rule','Resp'})
xlim([0.4 19.6])
set(gca,'FontSize',16,'LineWidth',2,'TickDir','out')
% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([analysis_figures_dir '\fMRI_decoding_',date,'.pdf'], '-dpdf', pdf_print_resolution)

