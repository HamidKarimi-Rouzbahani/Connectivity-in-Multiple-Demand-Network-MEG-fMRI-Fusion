clc
clear all;
close all;


Stim_Resp='Stim';
time_span=1; % time_span1: 2500;  time_span2: 3000;

analysis_figures_dir  = ['Z:\projects\Hamid\Projects\MD\Analyses\Results_temp_playing\Figures\'];

% if time_span==1
%     load(['MVGC_All_',Stim_Resp,'no_smoothing.mat'],'signss','pvals')
% else
%     load(['MVGC_All_',Stim_Resp,'no_smoothing.mat'],'signss','pvals')
% end
if time_span==1
    load(['MVGC_All_',Stim_Resp,'2500.mat'],'signss','pvals')
else
    load(['MVGC_All_',Stim_Resp,'3000.mat'],'signss','pvals')
end

Regions_summereized_reordered={'Visual','LOC','IPS','AI/FO','IFS','ACC'};
titles={'Coarse stimulus','Fine stimulus','Rule','Response'};
significance=0.05;


imagesc(squeeze(pvals(4,:,:)),[0 1])
%       imagesc(squeeze(signss(4,:,:)),[0 1])
%     axis equal
xticks([1:24])
xtickangle(45)
xticklabels(repmat(Regions_summereized_reordered,[1 4]))
xlabel('Source')
yticks([1:24])
ytickangle(45)
yticklabels(repmat(Regions_summereized_reordered,[1 4]))
ylabel('Destination')
CB=colorbar;
for i=1:size(pvals,2)
    for j=1:size(pvals,3)
        if ~isnan(pvals(cases,i,j))
            text(j-0.4,i,num2str(pvals(cases,i,j),'%1.1f'),'FontSize',8);
        end
    end
end
title(titles(cases))
%     yticks('')
%     xticks('')
set(gca,'fontsize', 10);

pdf_paper_size=[20 20];
pdf_print_resolution   = '-r300';
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
fig.PaperSize       = pdf_paper_size;
% print([analysis_figures_dir '\IFA_',Stim_Resp,'_',date,'.pdf'], '-dpdf', pdf_print_resolution)

