clear
close all
clc

param.data_path             = ['D:\Hamid\Postdoc\MD\Analyses\Results_temp_playing\'];
param.analysis_figures_dir  = ['D:\Hamid\Postdoc\MD\Analyses\Results_temp_playing\Figures\'];


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
        BFactors_stm(time_stim,iCat) = bf.ttest(param.aligned(1).data(:, time_stim, iCat),nanmean(param.aligned(1).data(:, 1:100, iCat),2));
    end
    for time_dec=1:size(param.aligned(2).data(:, :, iCat),2)
        BFactors_dec(time_dec,iCat) = bf.ttest(param.aligned(2).data(:, time_dec, iCat),nanmean(param.aligned(1).data(:, 1:100, iCat),2));
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
axis_ylim              = [45 75];
significant_points_y   = 48;
axis_ylim_Bayes         = [-1 13];
axis_ylim_spc          = 7;
axis_xlim_spc_Bayes       = 6;
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
pdf_file_name          = ['MEG_decoding_' date];
pdf_paper_size         = [20 20];
pdf_print_resolution   = '-r300';

%
% cl                     = colormap(hot);
% cl                     = cl(1:10:end, :);
cl=[[0.9 0.1 0.1];[0.1 0.9 0.1];[0.1 0.1 0.9];[0.5 0.5 0.5]];

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
span_stim=param.window_stim(1):param.slidwind:param.window_stim(2);
span_dec=param.window_dec(1):param.slidwind:param.window_dec(2);
ic = 1;
for iCat = [1 2 4 3]
    
    % for stim aligned
    this_data_stim = nanmean(param.aligned(1).data(:, :, iCat));
    this_data_stim_sd = nanstd(param.aligned(1).data(:, :, iCat));
    
    mean_data_stim = this_data_stim(:, param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2));
    mean_data_stim = smooth(mean_data_stim);

    
    this_data_dec = nanmean(param.aligned(2).data(:, :, iCat));
    this_data_dec_sd = nanstd(param.aligned(2).data(:, :, iCat));
    
    mean_data_dec = this_data_dec(:, param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2));
    mean_data_dec = smooth(mean_data_dec);

          
    BFactors_stm_tmp=log10(BFactors_stm(param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2),iCat));

    if ~isempty(find([smooth(BFactors_stm_tmp,Bayes_smoothing)>log10(param.up_thresh)]' & span_stim>30,1))
        tmp=find([smooth(BFactors_stm_tmp,Bayes_smoothing)>log10(param.up_thresh)]' & span_stim>30,1);
        first_significant_time(ic)=span_stim(tmp(1));
    else
        first_significant_time(ic)=nan;
    end
        
    for subject=1:24
        tmp_data=param.aligned(1).data(subject,param.time_stim >= param.window_stim(1) & param.time_stim <= param.window_stim(2), iCat);
        [~,max_time]=max(tmp_data(span_stim>30));
        peak_time_stim(subject,ic)=span_stim(max_time+(length(span_stim)-sum(span_stim>30)));
    end
    
    for subject=1:24
        tmp_data=param.aligned(2).data(subject,param.time_dec >= param.window_dec(1) & param.time_dec <= param.window_dec(2), iCat);
        [~,max_time]=max(tmp_data(span_dec<0));
        peak_time_dec(subject,ic)=abs(span_dec(max_time));
    end
    
    
    ic = ic + 1;        
end


    barwidths=0.3;
    colors_summerized={[0.1 0.7 0.1],[0.7 0.1 0.1],[0.1 0.1 0.7],[0.4 0.4 0.4]};


figure;
for iCat=1:4
    subplot(121)
    xtickpositions_smrzd=[1:barwidths*10:length(first_significant_time)*10*barwidths]+iCat*2*barwidths;
    bar(xtickpositions_smrzd,nanmean(peak_time_dec),'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
    hold on;
%     error(xtickpositions_smrzd,nanmean(peak_time_dec(:,iCat)),'facecolor',colors_summerized{iCat},'barwidth',barwidths*0.66,'EdgeColor','w');
end
% subplot(121)
% xticks([''])
% yticks([0:400:3000])
% yticklabels({'0','400','800','1200','1600','2000'});
% % ylim([0 2200])
% ylabel('Time to maximum decoding (ms)')
% box off;
% % xlim([0.4 19.6])
% xticks([1:barwidths*10:length(first_significant_time)*10*barwidths]+2.5*2*barwidths)
% xtickangle(45)
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','out')


% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
print([param.analysis_figures_dir '\MEG_decoding_stats2_',date,'.pdf'], '-dpdf', pdf_print_resolution)

