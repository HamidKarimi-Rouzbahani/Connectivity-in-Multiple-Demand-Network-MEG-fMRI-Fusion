clc;
close all;
clear all;
param.data_path             = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\'];
param.analysis_figures_dir  = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\Figures\'];

smoothing_rate_orig=50;
significance_level=0.001;



infos={'StimSide','StimPeriph','Rule','Response'};
Stim_Resp={'Stim','Resp'};
for OnsetResponse=1:2
    clear cors_summerizedNew
    for info=1:4
        load(['Fusion_',infos{info},'_',Stim_Resp{OnsetResponse},'_aligned.mat']);
        iteration=1;
        if smoothing_rate_orig>1
            for count=1:11
                for it=1:1000
                    fusion(:,count,it)=smooth(squeeze(fusion(:,count,it)),smoothing_rate_orig);
                end
            end
        end
        for region=1:6
            regions=[1 1;2:3; 4:5; 6:7; 8:9; 10:11];
            cors_summerizedNew(:,info,region)=nanmean(fusion(:,regions(region,:),iteration),2);
        end
    end
    param.aligned(OnsetResponse).data=cors_summerizedNew;
end
for OnsetResponse=1:2
    clear cors_summerizedNew significance
    if OnsetResponse==1
        significance=nan(900,4,6);
    else
        significance=nan(1000,4,6);
    end
    for info=1:4
        load(['Fusion_',infos{info},'_',Stim_Resp{OnsetResponse},'_aligned.mat']);
        if smoothing_rate_orig>1
            for count=1:11
                for it=1:1000
                    fusion(:,count,it)=smooth(squeeze(fusion(:,count,it)),smoothing_rate_orig);
                end
            end
        end
        iteration1=1;
        iteration2=2:1000;
        for region=1:6
            regions=[1 1;2:3; 4:5; 6:7; 8:9; 10:11];
            for time=1:length(nanmean(fusion(:,regions(region,:),iteration1),2))
                if nanmean(fusion(time,regions(region,:),iteration1),2)>0
                    if round((1-significance_level)*length(iteration2))<sum(nanmean(fusion(time,regions(region,:),iteration1),2)>nanmean(fusion(time,regions(region,:),iteration2),2))
                        significance(time,info,region)=1;
                    end
%                 else
%                     if round((1-significance_level)*length(iteration2))<sum(nanmean(fusion(time,regions(region,:),iteration1),2)<nanmean(fusion(time,regions(region,:),iteration2),2))
%                         significance(time,info,region)=1;
%                     end
                end
            end
            %% clustering significance
            spann=10;
            significance2=significance;
            for time=1:length(nanmean(fusion(:,regions(region,:),iteration1),2))-spann
                if nansum(significance(time:time+spann,info,region))<spann
                    significance2(time:time+spann-1,info,region)=nan;
                end
            end
            significance=significance2;
        end
    end
    param.aligned(OnsetResponse).significance=significance;
end


%% visualization
close all
clc
smoothing_rate_plot=1;

param.sr               = 1000; % sampling arte
param.window_stim      = [-200 2500]; % window of presentation
param.window_dec       = [-2500 500]; % window of presentation
param.slidwind         = 5;
param.time_stim        = [-500:param.slidwind :3995]+(smoothing_rate_orig*5)./2+10;
param.time_dec         = [-4000:param.slidwind :995]-(smoothing_rate_orig*5)./2+10;
param.subj_name        = 'all';
param.error            = 'sem';
param.up_thresh          = 3;
param.down_thresh          = param.up_thresh./10;
param.linestyle        = {'-','-','-', '-','-'};
% set plot properties
plot_linewidth         = 1;
plot_linewidth_added   =1.8;
plot_linewidth_Bayes         = 0.6;


% set axis properties
axis_ylim              = [-0.5 1.3];
axis_ylim_spc          = 7;
% axis_xlim              = param.window;
axis_tick_len          = 2;
axis_box_outline       = 'off';
axis_tick_style        = 'out';
axis_xlabel            = 'Time (ms)';
axis_ylabel            = {'Commonality index'};
axis_yxlabel_fontsize  = 10;
axis_xlabel_fontangle = 'normal';
axis_ylabel_fontangle = 'normal';
axis_yxtick_fontsize   = 8;
axis_yxtick_fontangle  = 'italic';
axis_title_fontangle   = 'normal';
axis_linewidth         = 1;

% set legend properties
legend_box_outline     = 'off';
legend_box_loaction    = 'northeast';
legend_fontangel       = 'normal';
legend_fontsize        = 8;

% set saveing/printing properties
pdf_file_name          = ['Fusion_' date];
pdf_paper_size         = [20 20];
pdf_print_resolution   = '-r300';

% cl                     = colormap(hot);
% cl                     = cl(1:10:end, :);
% cl=[[0.7 0.7 0.1];[0.7 0.1 0.7];[0.1 0.7 0.7];[0.1 0.7 0.7]];
% cl=[[0.7 0.1 0.1];[0.1 0.7 0.1];[0.1 0.1 0.7];[0.2 0.2 0.2]];
cl=colormap(hot(13));
col_order = [1:6];

gray_scale             = 0.7;

% categories:

iCat=1
Categories={'StimSide','StimPeriph','Rule','Response'};

Regions={'ACC','IPS','IFS','AI-FO','LOC','Visual'};
span_stim=param.window_stim(1):param.slidwind:param.window_stim(2);
span_dec=param.window_dec(1):param.slidwind:param.window_dec(2);
peak_time_dec=ones(6,4);
peak_time_stim=ones(6,4);
for iCat=1:4
    first_significant_time=ones(6,1);
    
    ic = 0;
    for region=[6 5 2 4 3 1]
        ic=ic+1;
        % for stim aligned
        this_data_stim = squeeze(param.aligned(1).data(:, iCat,region));
        this_data_stim = smooth(this_data_stim,smoothing_rate_plot);
        mean_data_stim = this_data_stim(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));

        this_significance_stim = squeeze(param.aligned(1).significance(:,iCat,region));
        significance_stim = this_significance_stim(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
        ind_first_sig_tmp=find(significance_stim(span_stim>30)==1);
        if ~isempty(ind_first_sig_tmp)
            ind_first_sig=ind_first_sig_tmp(1);
            first_sig_stim(ic,iCat)=span_stim(ind_first_sig+(length(span_stim)-sum(span_stim>30)));            
        else
            first_sig_stim(ic,iCat)=nan;           
        end
        
        [maxi,max_time]=max(mean_data_stim(span_stim>30));
        peak_time_stim(ic,iCat)=span_stim(max_time+(length(span_stim)-sum(span_stim>30)));
        maximum_stim(ic)=maxi;
               
        this_data_dec = squeeze(param.aligned(2).data(:, iCat,region));
        this_data_dec = smooth(this_data_dec,smoothing_rate_plot);
        mean_data_dec = this_data_dec(param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
        
%         this_significance_dec = squeeze(param.aligned(2).significance(:,iCat,region));
%         significance_dec = this_significance_dec(param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
%         ind_first_sig=find(significance_dec(span_dec<0));
%         ind_first_sig=ind_first_sig(1);
      
        
        [maxi,max_time]=max(mean_data_dec(span_dec<0));
        peak_time_dec(ic,iCat)=abs(span_dec(max_time));
        maximum_dec(ic)=maxi;
        
        this_data_stim_sig = squeeze(param.aligned(1).significance(:, iCat,region));
        data_stim_sig = this_data_stim_sig(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
        
        if ~isempty(find([data_stim_sig==1]' & span_stim>30,1))
            tmp=find([data_stim_sig==1]' & span_stim>30,1);
            first_significant_time(ic)=span_stim(tmp(1));
        else
            first_significant_time(ic)=nan;
        end
    end
    barwidths=0.3;
    colors_summerized={[0.7 0.1 0.1],[0.1 0.7 0.1],[0.1 0.1 0.7],[0.4 0.4 0.4]};
    xtickpositions_smrzd=[1:barwidths*10:length(first_significant_time)*10*barwidths]+iCat*2*barwidths;
    subplot(121)
    bar(xtickpositions_smrzd,maximum_stim,'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on;
    subplot(122)
    bar(xtickpositions_smrzd,maximum_dec,'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on;
end

% maximum commonality
subplot(121)
yticks(linspace(-0.05,1.3, 6))
% yticklabels({'0','200','400','600','800','1000'});
ylim([-0.06 1.3])
ylabel('Maximum commonality index')
box off;
xlim([0.4 19.6])
Regions={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
xticks([1:barwidths*10:length(first_significant_time)*10*barwidths]+2.5*2*barwidths)
xticklabels(Regions)
xtickangle(45)
plot([0.4 19.6],[0 0],'k','linewidth',1)
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','out')

subplot(122)
% yticks([0:400:2200])
yticks(linspace(-0.05,1.3, 6))
% yticklabels({'0','400','800','1200','1600','2000'});
% ylim([0 2200])
ylim([-0.06 1.3])
% ylabel('Maximum commonality index')
ylabel('Maximum commonality index')
box off;
xlim([0.4 19.6])
xticks([1:barwidths*10:length(first_significant_time)*10*barwidths]+2.5*2*barwidths)
Regions={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
xticklabels(Regions)
xtickangle(45)
plot([0.4 19.6],[0 0],'k','linewidth',1)
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','out')

% % set the prining properties
% fig                 = gcf;
% fig.PaperUnits      = 'centimeters';
% fig.Position        = [100 100 570 300];
% fig.PaperSize       = pdf_paper_size;
% % print([param.analysis_figures_dir '\Fusion_stats1_',date,'.pdf'], '-dpdf', pdf_print_resolution)
close all


% time to max
figure
% time to first significanace
% figure
subplot(131)
for iCat=1:4
    xtickpositions_smrzd=[1:barwidths*10:length(first_significant_time)*10*barwidths]+iCat*2*barwidths;
%     bar(xtickpositions_smrzd,peak_time_stim(:,iCat),'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    bar(xtickpositions_smrzd,first_sig_stim(:,iCat),'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on;hold on;
end
yticks([0:250:1000])
yticklabels({'0','250','750','1000'});
ylim([0 1000])
ylabel([{'Time to first'};{'significant commonality (ms)'}])
box off;
xlim([0.4 19.6])
xticks([1:barwidths*10:length(first_significant_time)*10*barwidths]+2.5*2*barwidths)
Regions={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
xticklabels(Regions)
xtickangle(45)
plot([0.4 19.6],[0 0],'k','linewidth',1)
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','out')


subplot(132)
for iCat=1:4
    xtickpositions_smrzd=[1:barwidths*10:length(first_significant_time)*10*barwidths]+iCat*2*barwidths;
    bar(xtickpositions_smrzd,peak_time_stim(:,iCat),'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on;
end
yticks([0:400:2200])
yticklabels({'0','400','800','1200','1600','2000'});
ylim([0 2200])
% ylabel('Time to maximum commonality (ms)')
ylabel([{'Time to'};{'maximum commonality (ms)'}])
box off;
xlim([0.4 19.6])
xticks([1:barwidths*10:length(first_significant_time)*10*barwidths]+2.5*2*barwidths)
Regions={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
xticklabels(Regions)
xtickangle(45)
plot([0.4 19.6],[0 0],'k','linewidth',1)
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','out')


subplot(133)
for iCat=1:4
    xtickpositions_smrzd=[1:barwidths*10:length(first_significant_time)*10*barwidths]+iCat*2*barwidths;
    bar(xtickpositions_smrzd,peak_time_dec(:,iCat),'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on;
end
yticks([0:400:2200])
yticklabels({'0','400','800','1200','1600','2000'});
ylim([0 2200])
% ylabel('Time since maximum commonality (ms)')
ylabel([{'Time since'};{'maximum commonality (ms)'}])

box off;
xlim([0.4 19.6])
xticks([1:barwidths*10:length(first_significant_time)*10*barwidths]+2.5*2*barwidths)
Regions={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
xticklabels(Regions)
xtickangle(45)
plot([0.4 19.6],[0 0],'k','linewidth',1)
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','out')


% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 50 570 300];
fig.PaperSize       = pdf_paper_size;
% print([param.analysis_figures_dir '\Fusion_statsAll_',date,'.pdf'], '-dpdf', pdf_print_resolution)

