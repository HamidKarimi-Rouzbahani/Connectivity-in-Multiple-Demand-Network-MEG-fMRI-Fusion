clear
close all
clc

param.data_path             = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\'];
param.analysis_figures_dir  = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\Figures\'];


infos={'Stim_side','Stim_periph','Rule','Resp'};
for OnsetResponse=1:2
    clear Accuracies
    for info=1:4
        for Subject=1:24
            if OnsetResponse==1
                load(['p',num2str(Subject),'_Stim_aligned_',infos{info},'_decoding_.mat'],'accuracy_correct');
%                 load(['p',num2str(Subject),'_Stim_aligned_',infos{info},'_decoding_svmlinear.mat'],'accuracy_correct');
            elseif OnsetResponse==2
                load(['p',num2str(Subject),'_Resp_aligned_',infos{info},'_decoding.mat'],'accuracy_correct');
%                 load(['p',num2str(Subject),'_Resp_aligned_',infos{info},'_decoding_svmlinear.mat'],'accuracy_correct');
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
categoriess=2; % 1=stim; 2=rule/resp

Bayes_smoothing   = 3;
panel_width=0.775;

param.window_stim      = [-500 2500]; % window of presentation
param.window_dec       = [-2500 500]; % window of presentation
param.slidwind         = 5;
param.time_stim        = -500:param.slidwind :3995;
param.time_dec         = -4000:param.slidwind :995;
param.up_thresh          = 3;
param.down_thresh          = 1/param.up_thresh;
param.linestyle        = {'-','-','-', '-','-'};



Sampling_rate_BF=1;
% set plot properties
plot_linewidth         = 0.7;
plot_linewidth_added   =1.3;
plot_linewidth_BF         = 1;


% set axis properties
if categoriess==1
    cats=[1 2];
    pdf_file_name          = ['MEG_decoding_Stims_' date];
    maxy_sig_BF=16;
    miny_sig_BF=-4;
    miny_sig_Decoding=45;
    maxy_sig_Decoding=75;
    x_lin_space_BF =12;
    y_lin_space_BF = 6;
    x_lin_space_Decoding = 12;
    y_lin_space_Decoding = 7;
else
    cats=[3 4];
    pdf_file_name          = ['MEG_decoding_RuleResp_' date];
    maxy_sig_BF=12;
    miny_sig_BF=-4;
    miny_sig_Decoding=47;
    maxy_sig_Decoding=62;
    x_lin_space_BF =12;
    y_lin_space_BF = 5;
    x_lin_space_Decoding = 12;
    y_lin_space_Decoding = 6;
end


axis_box_outline       = 'off';
xlabel_all            = 'Time (ms)';
ylabel_decoding           = 'Decoding accuracy (%)';
ylabel_BF            = 'BF (Log_{10})';
xylabel_fontsize  = 10;
axis_linewidth         = 1;


minx_plot_stim=-200;
maxx_plot_stim=2000;
minx_plot_resp=-1700;
maxx_plot_resp=500;

% set saveing/printing properties
pdf_paper_size         = [20 20];
pdf_print_resolution   = '-r300';

col_offset=0.3;
cl=[[0.9 col_offset col_offset];[col_offset 0.9 col_offset];[col_offset col_offset 0.9];[0.5 0.5 0.5];[0.9 0.9 col_offset]];

% Categories:
%     1- 'Stim side',
%     2- 'Stim periph',
%     3- 'Rule',
%     4- 'Resp',

for iCat = [1 2]
    this_data_stim = param.aligned(1).data(:, :, iCat);
%     this_data_stim = this_data_stim(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    
    this_data_dec = param.aligned(2).data(:, :, iCat);
%     this_data_dec = this_data_dec(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    
    data_stim(:,:,iCat)=this_data_stim;
    data_dec(:,:,iCat)=this_data_dec;
end
for time=1:size(data_stim,2)
    BF_coarse_fine_stim(time)=bf.ttest2(data_stim(:,time,1),data_stim(:,time,2));
end
for time=1:size(data_dec,2)
    BF_coarse_fine_dec(time)=bf.ttest2(data_dec(:,time,1),data_dec(:,time,2));
end

% extract the data
col_order = [1:5];

% Decoding stimulus-aligned
figure;
gca = axes('Position',[0.12 0.57 panel_width 0.4]); % xloc; ypos; width, height
for iCat = cats
    % for stim aligned
    this_data_stim = nanmean(param.aligned(1).data(:, :, iCat));
    this_data_stim_sd = nanstd(param.aligned(1).data(:, :, iCat));
    
    mean_data_stim = this_data_stim(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    mean_data_stim = smooth(mean_data_stim);
    
    std_data_stim = this_data_stim_sd(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    std_data_stim = smooth(std_data_stim);
    

    line([minx_plot_stim maxx_plot_stim],[50 50],'LineWidth',1,'Color','k','LineStyle',':');
    hold on;
    line([0 0],[miny_sig_Decoding maxy_sig_Decoding],'LineWidth',1,'Color','k','LineStyle',':');
    hh=shadedErrorBar(param.window_stim(1) : param.slidwind : param.window_stim(2),mean_data_stim,std_data_stim*1.96./sqrt(24),{'color',cl(col_order(iCat), :),'LineWidth',plot_linewidth},1);
    hh.LineWidth = plot_linewidth;
    hh.Color     = cl(col_order(iCat), :);
    hh.LineStyle = param.linestyle{1};
    hold on;
    BFactors_stm_tmp=smooth(log10(BFactors_stm(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2),iCat)),Bayes_smoothing);
    significant_points=nan(1,length(BFactors_stm_tmp));
    significant_points(BFactors_stm_tmp>log10(param.up_thresh))=1;
    hhh           = plot(param.window_stim(1) : param.slidwind : param.window_stim(2),mean_data_stim.*significant_points');
    hhh.LineWidth = plot_linewidth+plot_linewidth_added;
    hhh.Color     = cl(col_order(iCat), :);      
end
ylabel(ylabel_decoding)
box off;
xticklabels=num2str([linspace(minx_plot_stim,maxx_plot_stim,x_lin_space_Decoding)]');
xlabel(xlabel_all)

set(gca,'FontSize',xylabel_fontsize,'LineWidth',axis_linewidth,'XTick',...
    [linspace(minx_plot_stim,maxx_plot_stim,x_lin_space_Decoding)],'XTickLabel',...
    xticklabels,'YTick',...
    [linspace(miny_sig_Decoding,maxy_sig_Decoding,y_lin_space_Decoding)],'YTickLabel',{[linspace(miny_sig_Decoding,maxy_sig_Decoding,y_lin_space_Decoding)]},...
    'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k');
xlim([minx_plot_stim maxx_plot_stim])
ylim([miny_sig_Decoding maxy_sig_Decoding])
% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([ param.analysis_figures_dir '\' pdf_file_name '_stim.pdf'], '-dpdf', pdf_print_resolution)
    
% Decoding response-aligned

figure;
gca = axes('Position',[0.12 0.57 panel_width 0.4]); % xloc; ypos; width, height
for iCat=cats
    this_data_dec = nanmean(param.aligned(2).data(:, :, iCat));
    this_data_dec_sd = nanstd(param.aligned(2).data(:, :, iCat));
    
    mean_data_dec = this_data_dec(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    mean_data_dec = smooth(mean_data_dec);
    
    std_data_dec = this_data_dec_sd(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    std_data_dec = smooth(std_data_dec);
    
    line([minx_plot_resp maxx_plot_resp],[50 50],'LineWidth',1,'Color','k','LineStyle',':');
    hold on;
    line([0 0],[miny_sig_Decoding maxy_sig_Decoding],'LineWidth',1,'Color','k','LineStyle',':');
    hh=shadedErrorBar(param.window_dec(1) : param.slidwind : param.window_dec(2),mean_data_dec,std_data_dec*1.96./sqrt(24),{'color',cl(col_order(iCat), :),'LineWidth',plot_linewidth},1);
    hh.LineWidth = plot_linewidth;
    hh.Color     = cl(col_order(iCat), :);
    hh.LineStyle = param.linestyle{1};
    hold on;
    BFactors_dec_tmp=smooth(log10(BFactors_dec(param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2),iCat)),Bayes_smoothing);
    significant_points=nan(1,length(BFactors_dec_tmp));
    significant_points(BFactors_dec_tmp>log10(param.up_thresh))=1;
    hhh           = plot(param.window_dec(1) : param.slidwind : param.window_dec(2),mean_data_dec.*significant_points');
    hhh.LineWidth = plot_linewidth+plot_linewidth_added;
    hhh.Color     = cl(col_order(iCat), :);
end    
ylabel(ylabel_decoding)
box off;
xticklabels=num2str([linspace(minx_plot_resp,maxx_plot_resp,x_lin_space_Decoding)]');
xlabel(xlabel_all)

set(gca,'FontSize',xylabel_fontsize,'LineWidth',axis_linewidth,'XTick',...
    [linspace(minx_plot_resp,maxx_plot_resp,x_lin_space_Decoding)],'XTickLabel',...
    xticklabels,'YTick',...
    [linspace(miny_sig_Decoding,maxy_sig_Decoding,y_lin_space_Decoding)],'YTickLabel',{[linspace(miny_sig_Decoding,maxy_sig_Decoding,y_lin_space_Decoding)]},...
    'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k');
xlim([minx_plot_resp maxx_plot_resp])
ylim([miny_sig_Decoding maxy_sig_Decoding])

% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([ param.analysis_figures_dir '\' pdf_file_name '_resp.pdf'], '-dpdf', pdf_print_resolution)



% Bayes factor stimulus-aligned

if categoriess==1
    Effects=[BFactors_stm BF_coarse_fine_stim'];
    added_cat=5;
else
    Effects=[BFactors_stm];
    added_cat=[];
end
ii=1;
figure;
for iCat=[cats added_cat]   

    Null_color=[0 0 0];
    data_time_samples=[param.window_stim(1):5:param.window_stim(2)];
    Bayes_samples=1:Sampling_rate_BF:length(data_time_samples);  % revision 2

    gca = axes('Position',[0.12 0.8-(ii-1)*.22 panel_width 0.17]); % xloc; ypos; width, height 
    line([minx_plot_stim maxx_plot_stim],[0 0],'LineWidth',1,'Color','k','LineStyle',':');
    hold on;
    line([0 0],[miny_sig_BF maxy_sig_BF],'LineWidth',1,'Color','k','LineStyle',':');
    
    c=0;
    g=0;
    gg=0;
    ggg=0;
    for dots=data_time_samples
        c=c+1;
        if Effects(param.time_stim==dots,iCat)>=param.up_thresh
            plot(data_time_samples(c),log10(Effects(param.time_stim==dots,iCat)),'LineStyle','none','marker','o','MarkerFaceColor',cl(col_order(iCat),:),'Color',cl(col_order(iCat),:),'linewidth',.3,'markersize',3);
        elseif Effects(param.time_stim==dots,iCat)<param.up_thresh && Effects(param.time_stim==dots,iCat)>=param.down_thresh
            plot(data_time_samples(c),log10(Effects(param.time_stim==dots,iCat)),'LineStyle','none','marker','o','MarkerFaceColor',[1 1 1],'Color',Null_color,'linewidth',.3,'markersize',3);
        elseif Effects(param.time_stim==dots,iCat)<=param.down_thresh
            plot(data_time_samples(c),log10(Effects(param.time_stim==dots,iCat)),'LineStyle','none','marker','o','MarkerFaceColor',Null_color,'Color',Null_color,'linewidth',.3,'markersize',3);
        end
        hold on;
    end
    
    ylabel(ylabel_BF)
    box off;
    if (categoriess==1 && iCat==5) || (categoriess~=1 && iCat==4)
        xticklabels=num2str([linspace(minx_plot_stim,maxx_plot_stim,x_lin_space_BF)]');
        xlabel(xlabel_all)
    else
        xticklabels={};
    end
    
    set(gca,'FontSize',xylabel_fontsize,'LineWidth',axis_linewidth,'XTick',...
        [linspace(minx_plot_stim,maxx_plot_stim,x_lin_space_BF)],'XTickLabel',...
        xticklabels,'YTick',...
        [linspace(miny_sig_BF,maxy_sig_BF,y_lin_space_BF)],'YTickLabel',{[linspace(miny_sig_BF,maxy_sig_BF,y_lin_space_BF)]},...
        'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k');
    xlim([minx_plot_stim maxx_plot_stim])
    ylim([miny_sig_BF maxy_sig_BF])
    ii=ii+1;
end
% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([ param.analysis_figures_dir '\' pdf_file_name '_stim_BF.pdf'], '-dpdf', pdf_print_resolution)

% Bayes factor response-aligned
if categoriess==1
    Effects=[BFactors_dec BF_coarse_fine_dec'];
    added_cat=5;
else
    Effects=[BFactors_dec];
    added_cat=[];
end
ii=1;
figure;
for iCat=[cats added_cat]   

    data_time_samples=[param.window_dec(1):5:param.window_dec(2)];
    Bayes_samples=1:Sampling_rate_BF:length(data_time_samples);  

    gca = axes('Position',[0.12 0.8-(ii-1)*.22 panel_width 0.17]); % xloc; ypos; width, height 
    line([minx_plot_stim maxx_plot_stim],[0 0],'LineWidth',1,'Color','k','LineStyle',':');
    hold on;
    line([0 0],[miny_sig_BF maxy_sig_BF],'LineWidth',1,'Color','k','LineStyle',':');
    
    c=0;
    for dots=data_time_samples
        c=c+1;
        if Effects(param.time_dec==dots,iCat)>=param.up_thresh
            plot(data_time_samples(c),log10(Effects(param.time_dec==dots,iCat)),'LineStyle','none','marker','o','MarkerFaceColor',cl(col_order(iCat),:),'Color',cl(col_order(iCat),:),'linewidth',.3,'markersize',3);
        elseif Effects(param.time_dec==dots,iCat)<param.up_thresh && Effects(param.time_dec==dots,iCat)>=param.down_thresh
            plot(data_time_samples(c),log10(Effects(param.time_dec==dots,iCat)),'LineStyle','none','marker','o','MarkerFaceColor',[1 1 1],'Color',Null_color,'linewidth',.3,'markersize',3);
        elseif Effects(param.time_dec==dots,iCat)<param.down_thresh
            plot(data_time_samples(c),log10(Effects(param.time_dec==dots,iCat)),'LineStyle','none','marker','o','MarkerFaceColor',Null_color,'Color',Null_color,'linewidth',.3,'markersize',3);
        end
        hold on;
    end
    
    ylabel(ylabel_BF)
    box off;
    if (categoriess==1 && iCat==5) || (categoriess~=1 && iCat==4)
        xticklabels=num2str([linspace(minx_plot_resp,maxx_plot_resp,x_lin_space_BF)]');
        xlabel(xlabel_all)
    else
        xticklabels={};
    end
    
    set(gca,'FontSize',xylabel_fontsize,'LineWidth',axis_linewidth,'XTick',...
        [linspace(minx_plot_resp,maxx_plot_resp,x_lin_space_BF)],'XTickLabel',...
        xticklabels,'YTick',...
        [linspace(miny_sig_BF,maxy_sig_BF,y_lin_space_BF)],'YTickLabel',{[linspace(miny_sig_BF,maxy_sig_BF,y_lin_space_BF)]},...
        'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k');
    xlim([minx_plot_resp maxx_plot_resp])
    ylim([miny_sig_BF maxy_sig_BF])
    ii=ii+1;
end


% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([ param.analysis_figures_dir '\' pdf_file_name '_resp_BF.pdf'], '-dpdf', pdf_print_resolution)

