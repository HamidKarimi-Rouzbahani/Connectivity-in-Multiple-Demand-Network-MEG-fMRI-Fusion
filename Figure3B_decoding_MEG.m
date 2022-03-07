clear
close all
clc

param.data_path             = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\'];
param.analysis_figures_dir  = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\Figures\New\'];


infos={'Stim_side','Stim_periph','Rule','Resp'};
for OnsetResponse=1:2
    clear Accuracies
    for info=1:4
        for Subject=1:24
            if OnsetResponse==1
                load(['p',num2str(Subject),'_Stim_aligned_',infos{info},'_decoding_.mat'],'accuracy_correct');
            elseif OnsetResponse==2
                load(['p',num2str(Subject),'_Resp_aligned_',infos{info},'_decoding.mat'],'accuracy_correct');
            end
            Accuracies(Subject,:,info)=accuracy_correct*100;
        end
    end
    param.aligned(OnsetResponse).data=Accuracies;
end

for iCat = [1:4]
    for time_stim=1:size(param.aligned(1).data(:, :, iCat),2)
        BFactors_stm(time_stim,iCat) = bf.ttest2(param.aligned(1).data(:, time_stim, iCat),nanmean(param.aligned(1).data(:, 1:100, iCat),2));
    end
    for time_dec=1:size(param.aligned(2).data(:, :, iCat),2)
        BFactors_dec(time_dec,iCat) = bf.ttest2(param.aligned(2).data(:, time_dec, iCat),nanmean(param.aligned(1).data(:, 1:100, iCat),2));
    end
end

%% visualization
close all
clc
Bayes_smoothing   = 3;

param.sr               = 1000; % sampling arte
param.window_stim      = [-200 2500]; % window of presentation
param.window_dec       = [-2500 500]; % window of presentation
param.slidwind         = 5;
param.time_stim        = -500:param.slidwind :3995;
param.time_dec         = -4000:param.slidwind :995;
param.subj_name        = 'all';
param.error            = 'sem';
param.up_thresh          = 3;
param.down_thresh          = param.up_thresh./10;
param.linestyle        = {'-','-','-', '-','-'};
% set plot properties
plot_linewidth         = 0.7;
plot_linewidth_added   =1.3;
plot_linewidth_Bayes         = 1;


% set axis properties
% axis_xlim              = param.window;
axis_tick_len          = 2;
axis_box_outline       = 'off';
axis_tick_style        = 'out';
axis_xlabel            = 'Time (ms)';
axis_ylabel            = 'Decoding accuracy (%)';
axis_ylabel_Bayes            = ['Log_1_0','(Bayes factor)'];
axis_yxlabel_fontsize  = 10;
axis_yxlabel_fontangle = 'normal';
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
pdf_paper_size         = [20 20];
pdf_print_resolution   = '-r300';

%
% cl                     = colormap(hot);
% cl                     = cl(1:10:end, :);
cl=[[0.9 0.1 0.1];[0.1 0.9 0.1];[0.5 0.5 0.5];[0.1 0.1 0.9]];

gray_scale             = 0.7;

% categories:
%     1- 'Stim side',
%     2- 'Stim periph',
%     3- 'Rule',
%     4- 'Resp',

for iCat = [1 2]
    this_data_stim = param.aligned(1).data(:, :, iCat);    
    this_data_stim = this_data_stim(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));

    this_data_dec = param.aligned(2).data(:, :, iCat);    
    this_data_dec = this_data_dec(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    
    data_stim(:,:,iCat)=this_data_stim;
    data_dec(:,:,iCat)=this_data_dec;
end    
for time=1:size(data_stim,2)
    BF_stim(time)=bf.ttest2(data_stim(:,time,1),data_stim(:,time,2));
end
for time=1:size(data_dec,2)
    BF_dec(time)=bf.ttest2(data_dec(:,time,1),data_dec(:,time,2));
end
    
% extract the data
col_order = [1 2 4 3];
difference_color=[0.8 0.8 0.2];
ic = 1;

categoriess=1; % 1=stim; 2=rule/resp
if categoriess==1
    cats=[1 2];
    pdf_file_name          = ['MEG_decoding_Stims_' date];
    axis_ylim              = [45 75];
    significant_points_y   = 48;
    axis_ylim_Bayes         = [-1 17];
    axis_ylim_spc          = 7;
    axis_xlim_spc_Bayes       = 6;
else
    cats=[4 3];
    pdf_file_name          = ['MEG_decoding_RuleResp_' date];
    axis_ylim              = [47 62];
    significant_points_y   = 48;
    axis_ylim_Bayes         = [-1 9.5];
    axis_ylim_spc          = 6;
    axis_xlim_spc_Bayes       = 6;
end

spann_clustering=0;

for iCat = cats
    ic=iCat;
    % for stim aligned
    this_data_stim = nanmean(param.aligned(1).data(:, :, iCat));
    this_data_stim_sd = nanstd(param.aligned(1).data(:, :, iCat));
    
    mean_data_stim = this_data_stim(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    mean_data_stim = smooth(mean_data_stim);
    
    std_data_stim = this_data_stim_sd(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    std_data_stim = smooth(std_data_stim);
    
    this_data_dec = nanmean(param.aligned(2).data(:, :, iCat));
    this_data_dec_sd = nanstd(param.aligned(2).data(:, :, iCat));
    
    mean_data_dec = this_data_dec(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    mean_data_dec = smooth(mean_data_dec);
    
    std_data_dec = this_data_dec_sd(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    std_data_dec = smooth(std_data_dec);
    
    
    subplot(2, 2, 3)
    if iCat<3
        ht           = plot(param.window_stim(1) : param.slidwind : param.window_stim(2),smooth(log10(BF_stim),Bayes_smoothing));
        ht.LineWidth = plot_linewidth_Bayes;
        ht.Color     = [0.8 0.8 0.1];
        ht.LineStyle = param.linestyle{1};
        hold on;
    end
    
    BFactors_stm_tmp=smooth(log10(BFactors_stm(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2),iCat)),Bayes_smoothing);
    h           = plot(param.window_stim(1) : param.slidwind : param.window_stim(2),BFactors_stm_tmp);
    h.LineWidth = plot_linewidth_Bayes;
    h.Color     = cl(col_order(ic), :);
    h.LineStyle = param.linestyle{1};
    hold on
    
    subplot(2, 2, 4)
    if iCat<3
        ht           = plot(param.window_dec(1) : param.slidwind : param.window_dec(2),smooth(log10(BF_dec),Bayes_smoothing));
        ht.LineWidth = plot_linewidth_Bayes;
        ht.Color     = difference_color;
        ht.LineStyle = param.linestyle{1};
        hold on;
    end
    
    BFactors_dec_tmp=smooth(log10(BFactors_dec(param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2),iCat)),Bayes_smoothing);
    h           = plot(param.window_dec(1) : param.slidwind : param.window_dec(2),BFactors_dec_tmp);
    h.LineWidth = plot_linewidth_Bayes;
    h.Color     = cl(col_order(ic), :);
    h.LineStyle = param.linestyle{1};
    hold on
    
    
    
    subplot(2, 2, 1)
    if iCat==1
        sig_difference=find(BF_stim>param.up_thresh);
        % clustering
        sig_difference2=nan(length(BF_stim),1);
        sig_difference2(sig_difference)=1;
        for time=1:length(sig_difference2)-spann_clustering
            if nansum(sig_difference2(time:time+spann_clustering-1))<spann_clustering
                sig_difference2(time)=nan;
            end
        end
        sig_difference=find(~isnan(sig_difference2));
        timings=param.window_stim(1) : param.slidwind : param.window_stim(2);
        plot(timings(sig_difference),repmat(significant_points_y,[1 length(timings(sig_difference))]),'.','color',difference_color);
        hold on;
    end
    % plot ocip
    hh=shadedErrorBar(param.window_stim(1) : param.slidwind : param.window_stim(2),mean_data_stim,std_data_stim*1.96./sqrt(24),{'color',cl(col_order(ic), :),'LineWidth',plot_linewidth},1);
    hh.LineWidth = plot_linewidth;
    hh.Color     = cl(col_order(ic), :);
    hh.LineStyle = param.linestyle{1};
    hold on;
    significant_points=nan(1,length(BFactors_stm_tmp));
    significant_points(BFactors_stm_tmp>log10(param.up_thresh))=1;
    hhh           = plot(param.window_stim(1) : param.slidwind : param.window_stim(2),mean_data_stim.*significant_points');
    hhh.LineWidth = plot_linewidth+plot_linewidth_added;
    hhh.Color     = cl(col_order(ic), :);
    hold on
    

    subplot(2, 2, 2)
    if iCat==1
        sig_difference=find(BF_dec>param.up_thresh);
        % clustering
        sig_difference2=nan(length(BF_stim),1);
        sig_difference2(sig_difference)=1;
        for time=1:length(sig_difference2)-spann_clustering
            if nansum(sig_difference2(time:time+spann_clustering-1))<spann_clustering
                sig_difference2(time)=nan;
            end
        end
        sig_difference=find(~isnan(sig_difference2));
        
        timings=param.window_dec(1) : param.slidwind : param.window_dec(2);
        plot(timings(sig_difference),repmat(significant_points_y,[1 length(timings(sig_difference))]),'.','color',difference_color);
        hold on;
    end
    hh=shadedErrorBar(param.window_dec(1) : param.slidwind : param.window_dec(2),mean_data_dec,std_data_dec*1.96./sqrt(24),{'color',cl(col_order(ic), :),'LineWidth',plot_linewidth},1);
    hh.LineWidth = plot_linewidth;
    hh.Color     = cl(col_order(ic), :);
    hh.LineStyle = param.linestyle{1};
    hold on;
    significant_points=nan(1,length(BFactors_dec_tmp));
    significant_points(BFactors_dec_tmp>log10(param.up_thresh))=1;           
    hhh           = plot(param.window_dec(1) : param.slidwind : param.window_dec(2),mean_data_dec.*significant_points');
    hhh.LineWidth = plot_linewidth+plot_linewidth_added;
    hhh.Color     = cl(col_order(ic), :);
    hold on
    
%     ic = ic + 1;
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
aX.XLabel.FontAngle = axis_yxlabel_fontangle;
aX.YLabel.String    = axis_ylabel;
aX.YLabel.FontSize  = axis_yxlabel_fontsize;
aX.YLabel.FontAngle = axis_yxlabel_fontangle;
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
aX.XLabel.FontAngle = axis_yxlabel_fontangle;
aX.YLabel.String    = axis_ylabel;
aX.YLabel.FontSize  = axis_yxlabel_fontsize;
aX.YLabel.FontAngle = axis_yxlabel_fontangle;
aX.YLim             = axis_ylim;
aX.YTick            = linspace(axis_ylim(1), axis_ylim(2), axis_ylim_spc);
aX.XTick            = linspace(param.window_dec(1), param.window_dec(end), 7);
aX.XTickLabel       = linspace(param.window_dec(1), param.window_dec(end), 7);
aX.YAxis.Visible    = 'off';
aX.LineWidth        = axis_linewidth;
aX.XLim             = [param.window_dec(1), param.window_dec(end)];

subplot(2, 2, 3)
plot(param.window_stim, log10(param.up_thresh)*[1 1], ':k')
plot(param.window_stim, log10(param.down_thresh)*[1 1], ':k')
plot([0 0], axis_ylim_Bayes , ':k'), hold on
% change the axis properties
aX                  = gca;
aX.Box              = axis_box_outline;
aX.TickDir          = axis_tick_style;
aX.TickLength       = axis_tick_len * aX.TickLength;
aX.FontSize         = axis_yxtick_fontsize;
aX.FontAngle        = axis_yxtick_fontangle;
aX.XLabel.String    = axis_xlabel;
aX.XLabel.FontSize  = axis_yxlabel_fontsize;
aX.XLabel.FontAngle = axis_yxlabel_fontangle;
aX.YLabel.String    = axis_ylabel_Bayes;
aX.YLabel.FontSize  = axis_yxlabel_fontsize;
aX.YLabel.FontAngle = axis_yxlabel_fontangle;
aX.YLim             = axis_ylim_Bayes;
aX.YTick            = linspace(axis_ylim_Bayes(1), axis_ylim_Bayes(2), axis_xlim_spc_Bayes);
aX.XTick            = linspace(param.window_stim(1), param.window_stim(end), 7);
aX.XTickLabel       = linspace(param.window_stim(1), param.window_stim(end), 7);
aX.LineWidth        = axis_linewidth;
aX.XLim             = [param.window_stim(1), param.window_stim(end)];

subplot(2, 2, 4)
plot(param.window_dec, log10(param.up_thresh)*[1 1], ':k')
plot(param.window_dec, log10(param.down_thresh)*[1 1], ':k')
plot([0 0], axis_ylim_Bayes , ':k'), hold on
% change the axis properties
aX                  = gca;
aX.Box              = axis_box_outline;
aX.TickDir          = axis_tick_style;
aX.TickLength       = axis_tick_len * aX.TickLength;
aX.FontSize         = axis_yxtick_fontsize;
aX.FontAngle        = axis_yxtick_fontangle;
aX.XLabel.String    = axis_xlabel;
aX.XLabel.FontSize  = axis_yxlabel_fontsize;
aX.XLabel.FontAngle = axis_yxlabel_fontangle;
% aX.YLabel.String    = axis_ylabel_Bayes;
% aX.YLabel.FontSize  = axis_yxlabel_fontsize;
% aX.YLabel.FontAngle = axis_yxlabel_fontangle;
aX.YAxis.Visible    = 'off';
aX.YLim             = axis_ylim_Bayes;
aX.YTick            = linspace(axis_ylim_Bayes(1), axis_ylim_Bayes(2), axis_xlim_spc_Bayes);
aX.XTick            = linspace(param.window_dec(1), param.window_dec(end), 7);
aX.XTickLabel       = linspace(param.window_dec(1), param.window_dec(end), 7);
aX.LineWidth        = axis_linewidth;
aX.XLim             = [param.window_dec(1), param.window_dec(end)];


% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([ param.analysis_figures_dir '\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)

