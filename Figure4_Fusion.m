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
param.up_thresh        = 3;
param.down_thresh          = param.up_thresh./10;
param.linestyle        = {'-','-','-', '-','-'};
% set plot properties
plot_linewidth         = 1;
plot_linewidth_added   =1.8;
plot_linewidth_Bayes   = 0.6;


% set axis properties
axis_ylim              = [-0.5 1.5];
axis_ylim_spc          = 5;
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


cl=[[0.7 0.1 0.1];[0.1 0.7 0.1];[0.1 0.1 0.7];[0.4 0.4 0.4]];

gray_scale             = 0.7;

% categories:

col_order = [1 2 4 3];
ic = 1;
region=4
Regions={'pre-SMA','IPS','IFS','AI-FO','LOC','Visual'};
for iCat = [1 2 4 3]
    
    % for stim aligned
    this_data_stim = squeeze(param.aligned(1).data(:, iCat,region));
    this_data_stim = smooth(this_data_stim,smoothing_rate_plot);
    mean_data_stim = this_data_stim(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    
    this_data_stim_sig = squeeze(param.aligned(1).significance(:, iCat,region));
    data_stim_sig = this_data_stim_sig(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    
    this_data_dec = squeeze(param.aligned(2).data(:, iCat,region));
    this_data_dec = smooth(this_data_dec,smoothing_rate_plot);
    mean_data_dec = this_data_dec(param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    
    this_data_dec_sig = squeeze(param.aligned(2).significance(:, iCat,region));
    data_dec_sig = this_data_dec_sig(param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    
    subplot(2, 2, 1)
    hh=plot(param.window_stim(1) : param.slidwind : param.window_stim(2),mean_data_stim);
    hh.LineWidth = plot_linewidth;
    hh.Color     = cl(col_order(ic), :);
    hh.LineStyle = param.linestyle{1};
    hold on;
    hhh           = plot(param.window_stim(1) : param.slidwind : param.window_stim(2),mean_data_stim.*data_stim_sig);
    hhh.LineWidth = plot_linewidth+plot_linewidth_added;
    hhh.Color     = cl(col_order(ic), :);
    hold on
    
    subplot(2, 2, 2)
    hh=plot(param.window_dec(1) : param.slidwind : param.window_dec(2),mean_data_dec);
    hh.LineWidth = plot_linewidth;
    hh.Color     = cl(col_order(ic), :);
    hh.LineStyle = param.linestyle{1};
    hold on;
    hhh           = plot(param.window_dec(1) : param.slidwind : param.window_dec(2),mean_data_dec.*data_dec_sig);
    hhh.LineWidth = plot_linewidth+plot_linewidth_added;
    hhh.Color     = cl(col_order(ic), :);
    hold on
    ic = ic + 1;
end

subplot(2, 2, 1)
plot(param.window_stim, 50*[1 1], ':k')
plot([0 0], axis_ylim , ':k'), hold on
plot([param.window_stim(1), param.window_stim(end)], [0 0], ':k')
% change the axis properties
aX                  = gca;
aX.Box              = axis_box_outline;
aX.TickDir          = axis_tick_style;
aX.TickLength       = axis_tick_len * aX.TickLength;
aX.FontSize         = axis_yxtick_fontsize;
aX.FontAngle        = axis_yxtick_fontangle;
aX.XLabel.String    = axis_xlabel;
aX.XLabel.FontSize  = axis_yxlabel_fontsize;
aX.XLabel.FontAngle = axis_xlabel_fontangle;
aX.YLabel.String    = axis_ylabel;
aX.YLabel.FontSize  = axis_yxlabel_fontsize;
aX.YLabel.FontAngle = axis_ylabel_fontangle;
aX.YLim             = axis_ylim;
aX.YTick            = linspace(axis_ylim(1), axis_ylim(2), axis_ylim_spc);
aX.XTick            = linspace(param.window_stim(1), param.window_stim(end), 7);
aX.XTickLabel       = linspace(param.window_stim(1), param.window_stim(end), 7);
aX.LineWidth        = axis_linewidth;
aX.XLim             = [param.window_stim(1), param.window_stim(end)];

subplot(2, 2, 2)
plot(param.window_dec, 50*[1 1], ':k')
plot([0 0], axis_ylim , ':k'), hold on
plot([param.window_dec(1), param.window_dec(end)], [0 0], ':k')
% change the axis properties
aX                  = gca;
aX.Box              = axis_box_outline;
aX.TickDir          = axis_tick_style;
aX.TickLength       = axis_tick_len * aX.TickLength;
aX.FontSize         = axis_yxtick_fontsize;
aX.FontAngle        = axis_yxtick_fontangle;
aX.XLabel.String    = axis_xlabel;
aX.XLabel.FontSize  = axis_yxlabel_fontsize;
aX.XLabel.FontAngle = axis_xlabel_fontangle;
aX.YLabel.String    = axis_ylabel;
aX.YLabel.FontSize  = axis_yxlabel_fontsize;
aX.YLabel.FontAngle = axis_ylabel_fontangle;
aX.YLim             = axis_ylim;
aX.YTick            = linspace(axis_ylim(1), axis_ylim(2), axis_ylim_spc);
aX.XTick            = linspace(param.window_dec(1), param.window_dec(end), 7);
aX.XTickLabel       = linspace(param.window_dec(1), param.window_dec(end), 7);
% aX.YAxis.Visible    = 'off';
aX.LineWidth        = axis_linewidth;
aX.XLim             = [param.window_dec(1), param.window_dec(end)];


% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([ param.analysis_figures_dir '\' ['Region_' Regions{region} pdf_file_name '.pdf']], '-dpdf', pdf_print_resolution)

% pause(10)
% end