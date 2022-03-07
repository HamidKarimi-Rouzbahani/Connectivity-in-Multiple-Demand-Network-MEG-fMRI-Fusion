%% Plotting across subejcts
clc;
clear all;
close all;
load('behavioural_data_MEG.mat')
load('behavioural_data_fMRI.mat')
upper_threshold=3;
lower_threshold=upper_threshold./10;

analysis_figures_dir  = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\Figures\New\'];

legend_box_outline     = 'off';
legend_box_loaction    = 'northeast';
legend_fontangel       = 'normal';
legend_fontsize        = 14;
pdf_paper_size=[20 20];
pdf_print_resolution   = '-r300';

widths=0.5;


 % trial_outcome: 1=stim error/coarse stim error; 2=both error/fine stim error; 3=rule error; 4= correct trials
colors_fMRI={[0.1 0.5 0.1],[0.5 0.1 0.1],[0.1 0.1 0.5],[0.3 0.3 0.3]};
colors_MEG={[0.1 0.9 0.1],[0.9 0.1 0.1],[0.1 0.1 0.9],[0.5 0.5 0.5]};

steps=0;
counter=0;
data=nan(2,4,30);
figure;
for trial_outcome=[2 1 3]
    steps=steps+2.5;
    counter=counter+1;

    subplot(121)
    plots(counter)=bar(steps,100-nanmean(Stims_subj(:,trial_outcome))*100,'facecolor',colors_fMRI{trial_outcome},'linewidth',widths,'EdgeColor','w');
    hold on;
    errorbar(steps,100-nanmean(Stims_subj(:,trial_outcome))*100,nanstd(Stims_subj(:,trial_outcome)*100)./(sqrt(30)/1.96),'linewidth',2,'color','k','capsize',0,'linestyle','none')
    
    
    plots(counter+3)=bar(steps+widths*1.8,nanmean(Acc_subj_MEG(:,trial_outcome))*100,'facecolor',colors_MEG{trial_outcome},'linewidth',widths,'EdgeColor','w');
    hold on;
    errorbar(steps+widths*1.8,nanmean(Acc_subj_MEG(:,trial_outcome))*100,nanstd(Acc_subj_MEG(:,trial_outcome)*100)./(sqrt(24)/1.96),'linewidth',2,'color','k','capsize',0,'linestyle','none')

    data(1,counter,1:30)=Stims_subj(:,trial_outcome);
    data(2,counter,1:24)=Acc_subj_MEG(:,trial_outcome);
end

combinations=[1 2;1 3;2 3];
for combs=1:3
    Bayesfactors_fMRI(combs)=bf.ttest(squeeze(data(1,combinations(combs,1),:)),squeeze(data(1,combinations(combs,2),:)));
    Bayesfactors_MEG(combs)=bf.ttest(squeeze(data(2,combinations(combs,1),:)),squeeze(data(2,combinations(combs,2),:)));
end
num_to_plot=round(Bayesfactors_fMRI(1),2,'significant');
if num_to_plot>upper_threshold || num_to_plot<lower_threshold
    text(2.5,100,num2str(num_to_plot,'%10.1e\n'),'color','k','FontWeight','bold','fontangle','italic');
else
    text(2.5,100,num2str(num_to_plot,'%10.1e\n'),'color','k','FontWeight','light','fontangle','normal');
end
num_to_plot=round(Bayesfactors_MEG(1),2,'significant');
if num_to_plot>upper_threshold || num_to_plot<lower_threshold
    text(3.5,101,num2str(num_to_plot,'%10.1e\n'),'color','k','FontWeight','bold','fontangle','italic');
else
    text(3.5,101,num2str(num_to_plot,'%10.1e\n'),'color','k','FontWeight','light','fontangle','normal');
end

subplot(121)
box off;
xticks([1:3])
ylabel('Behavioural accuracy (%)')
ylim([85 102])
xticks(['','',''])
set(gca,'FontSize',16,'LineWidth',2,'TickDir','out')

% % % % change legend properties
% h1=legend([plots(1),plots(4),plots(2),plots(5),plots(3),plots(6)],...
%     {'Stim Side (fMRI)','Stim Side (MEG)',...
%     'Stim Periphery (fMRI)','Stim Periphery (MEG)',...
%     'Rule (fMRI)','Rule (MEG)'});
% 
% h1.Location  = legend_box_loaction;
% h1.Box       = legend_box_outline;
% h1.FontAngle = legend_fontangel;
% h1.FontSize  = legend_fontsize;
% 



subplot(122)
trial_outcome=4;
plotss(1)=bar(1,nanmean(Stims_subj_rt(:,trial_outcome))*1000,'facecolor',colors_fMRI{trial_outcome},'linewidth',widths,'EdgeColor','w');
hold on;
errorbar(1,nanmean(Stims_subj_rt(:,trial_outcome))*1000,nanstd(Stims_subj_rt(:,trial_outcome)*1000)./(sqrt(30)/1.96),'linewidth',2,'color','k','capsize',0,'linestyle','none')

plotss(2)=bar(1+widths*1.8,nanmean(RT_subj_MEG(:,trial_outcome))*1000,'facecolor',colors_MEG{trial_outcome},'linewidth',widths,'EdgeColor','w');
hold on;
errorbar(1+widths*1.8,nanmean(RT_subj_MEG(:,trial_outcome))*1000,nanstd(RT_subj_MEG(:,trial_outcome)*1000)./(sqrt(30)/1.96),'linewidth',2,'color','k','capsize',0,'linestyle','none')
yticks([1250:250:2000])
box off;
xticks([''])
ylabel('Reaction time (ms)')
ylim([1250 2000])



% % % % change legend properties
% h2=legend([plotss(1),plotss(2)],{'fMRI','MEG'});
% h2.Location  = legend_box_loaction;
% h2.Box       = legend_box_outline;
% h2.FontAngle = legend_fontangel;
% h2.FontSize  = legend_fontsize;

set(gca,'FontSize',16,'LineWidth',2,'TickDir','out')
% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([analysis_figures_dir '\Behavioural_',date,'.pdf'], '-dpdf', pdf_print_resolution)
