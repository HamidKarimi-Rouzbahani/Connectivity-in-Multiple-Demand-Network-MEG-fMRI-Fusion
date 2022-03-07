clc;
clear all;
close all;

Behavioural_directory='D:\Hamid\Postdoc\MD\fMRI_Behavioural_data\';

Stims=nan*ones(4,4,6,31);
Stims_no=nan*ones(4,4,6,31);
Stims_rt=nan*ones(4,4,6,31);
Rules=nan*ones(4,4,6,31);
Rules_no=nan*ones(4,4,6,31);
Rules_rt=nan*ones(4,4,6,31);

s=0;
for subject=[1:26 28:31]
    runss=dir([Behavioural_directory,'/test_res_p',num2str(subject),'_r*.mat']);
    runs=[1:size(runss,1)];
    rule_version=[1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2];
    
    s=s+1;
    for run=runs
        clearvars -except runs s Stims Stims_no Stims_rt Rules Rules_no Rules_rt csub Behavioural_directory collapsed rule_version Subjects subject Main_analysis_directory run performances
        %     load([Directory,Subjects{subject},'/test_res_p',...
        %         num2str(subject),'_r',num2str(run),'.mat']);
        csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
        csub = csub(1:4);
        load([Behavioural_directory,'/test_res_p',num2str(subject),'_r',num2str(run),'.mat']);
        if subject==1
            first_stim_delay=4-test_res.allresults(1,12);
        else
            first_stim_delay=14-test_res.allresults(1,12);
        end
        all_responded_trials=find(~isnan(test_res.allresults(:,15)))';
        response_trial_outcome=zeros(length(test_res.allresults(:,12)),4);
        % Trial trial_outcome categorization
        tmp_ro=zeros(length(test_res.allresults(:,12)),4);
        
        for i=all_responded_trials
            
            if test_res.allresults(i,4)==rule_version(subject)
                
                if test_res.allresults(i,7)==1 && test_res.allresults(i,9)==1
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==2
                    % Stim + Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==3
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==4
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                    
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==1
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==2
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==3
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==4
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                    
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==1
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==2
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==3
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==4
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                    
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==1
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==2
                    % Corrrect
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==3
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==4
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                end
                
            else
                if test_res.allresults(i,7)==1 && test_res.allresults(i,9)==1
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==2
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==3
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==4
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                    
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==1
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==2
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==3
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==4
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                    
                    
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==1
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==2
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==3
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==4
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                    
                    
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==1
                    % Correct
                    tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,4)=i;
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==2
                    % Rule
                    tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,3)=i;
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==3
                    % Stim+Rule
                    tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,2)=i;
                elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==4
                    % Stim
                    tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                    response_trial_outcome(i,1)=i;
                end
            end
            %         Responded_RT_tmp(i,find(tmp_ro(i,:)~=0))=test_res.allresults(i,15);
        end
        
        for blk=1:2
            response_trial_outcome_blk=response_trial_outcome([[1:80]+(blk-1)*80],:);
            
            for trial_trial_outcome=1:4
                % Row1: Stim error;
                % Row2: Stim and Rule error;
                % Row3: Rule error;
                % Row4: Correct;
                Trials_outcome{trial_trial_outcome}=response_trial_outcome_blk(response_trial_outcome_blk(:,trial_trial_outcome)~=0,trial_trial_outcome);
                target_trials=Trials_outcome{trial_trial_outcome}';
                
                
                tmp_stim=zeros(size(response_trial_outcome_blk,1),4);
                tmp_rule=zeros(size(response_trial_outcome_blk,1),2);
                tmp_rule_color=zeros(size(response_trial_outcome_blk,1),4);
                tmp_resp=zeros(size(response_trial_outcome_blk,1),4);
                tmp_stim_rt=zeros(size(response_trial_outcome_blk,1),4);
                tmp_rule_rt=zeros(size(response_trial_outcome_blk,1),2);
                tmp_rule_color_rt=zeros(size(response_trial_outcome_blk,1),4);
                tmp_resp_rt=zeros(size(response_trial_outcome_blk,1),4);
                tmp_stim_col_16=zeros(size(response_trial_outcome_blk,1),16);
                tmp_stim_col_rt_16=zeros(size(response_trial_outcome_blk,1),16);
                tmp_stim_col_8=zeros(size(response_trial_outcome_blk,1),8);
                tmp_stim_col_rt_8=zeros(size(response_trial_outcome_blk,1),8);
                tmp_stim_col_4=zeros(size(response_trial_outcome_blk,1),4);
                tmp_stim_col_rt_4=zeros(size(response_trial_outcome_blk,1),4);
                
                for i=target_trials
                    t=i-(blk-1)*80;
                    tmp_stim(t,test_res.allresults(i,7))=test_res.allresults(i,12)+first_stim_delay;
                    tmp_stim_rt(t,test_res.allresults(i,7))=test_res.allresults(i,15);
                    
                    bias=(test_res.allresults(i,4)-1)*8+(test_res.allresults(i,5)-1)*4;
                    tmp_stim_col_16(t,test_res.allresults(i,7)+bias)=test_res.allresults(i,12)+first_stim_delay;
                    tmp_stim_col_rt_16(t,test_res.allresults(i,7)+bias)=test_res.allresults(i,15);
                    
                    bias=(test_res.allresults(i,4)-1)*4;
                    tmp_stim_col_8(t,test_res.allresults(i,7)+bias)=test_res.allresults(i,12)+first_stim_delay;
                    tmp_stim_col_rt_8(t,test_res.allresults(i,7)+bias)=test_res.allresults(i,15);
                    
                    bias=(test_res.allresults(i,4)-1)*2;
                    if test_res.allresults(i,7)>2
                        tmp=test_res.allresults(i,7)-2;
                    else
                        tmp=test_res.allresults(i,7);
                    end
                    tmp_stim_col_4(t,tmp+bias)=test_res.allresults(i,12)+first_stim_delay;
                    tmp_stim_col_rt_4(t,tmp+bias)=test_res.allresults(i,15);
                    
                    tmp_rule(t,test_res.allresults(i,4))=test_res.allresults(i,12)+first_stim_delay;
                    tmp_rule_rt(t,test_res.allresults(i,4))=test_res.allresults(i,15);
                    
                    if test_res.allresults(i,4)==1
                        tmp_rule_color(t,test_res.allresults(i,5))=test_res.allresults(i,12)+first_stim_delay;
                        tmp_rule_color_rt(t,test_res.allresults(i,5))=test_res.allresults(i,15);
                    else
                        tmp_rule_color(t,test_res.allresults(i,5)+test_res.allresults(i,4))=test_res.allresults(i,12)+first_stim_delay;
                        tmp_rule_color_rt(t,test_res.allresults(i,5)+test_res.allresults(i,4))=test_res.allresults(i,15);
                    end
                    
                    tmp_resp(t,test_res.allresults(i,9))=test_res.allresults(i,12)+first_stim_delay;
                    tmp_resp_rt(t,test_res.allresults(i,9))=test_res.allresults(i,15);
                end
                
                for cond=[1:4]
                    cond_shift=cond+(blk-1)*4;
                    Stim{trial_trial_outcome,cond_shift}=tmp_stim(tmp_stim(:,cond)~=0,cond);
                    Stim_rt{trial_trial_outcome,cond_shift}=tmp_stim_rt(tmp_stim_rt(:,cond)~=0,cond);
                    
                    if cond<3
                        cond_shiftR=cond+(blk-1)*2;
                        Rule{trial_trial_outcome,cond_shiftR}=tmp_rule(tmp_rule(:,cond)~=0,cond);
                        Rule_rt{trial_trial_outcome,cond_shiftR}=tmp_rule_rt(tmp_rule_rt(:,cond)~=0,cond);
                    end
                    Rule_color{trial_trial_outcome,cond_shift}=tmp_rule_color(tmp_rule_color(:,cond)~=0,cond);
                    Rule_color_rt{trial_trial_outcome,cond_shift}=tmp_rule_color_rt(tmp_rule_color_rt(:,cond)~=0,cond);
                    
                    Resp{trial_trial_outcome,cond_shift}=tmp_resp(tmp_resp(:,cond)~=0,cond);
                    Resp_rt{trial_trial_outcome,cond_shift}=tmp_resp_rt(tmp_resp_rt(:,cond)~=0,cond);
                end
                
                for cond=[1:16]
                    cond_shift=cond+(blk-1)*16;
                    Stim_Col_16{trial_trial_outcome,cond_shift}=tmp_stim_col_16(tmp_stim_col_16(:,cond)~=0,cond);
                    Stim_Col_rt_16{trial_trial_outcome,cond_shift}=tmp_stim_col_rt_16(tmp_stim_col_rt_16(:,cond)~=0,cond);
                end
                
                for cond=[1:8]
                    cond_shift=cond+(blk-1)*8;
                    Stim_Col_8{trial_trial_outcome,cond_shift}=tmp_stim_col_8(tmp_stim_col_8(:,cond)~=0,cond);
                    Stim_Col_rt_8{trial_trial_outcome,cond_shift}=tmp_stim_col_rt_8(tmp_stim_col_rt_8(:,cond)~=0,cond);
                end
                
                for cond=[1:4]
                    cond_shift=cond+(blk-1)*4;
                    Stim_Col_4{trial_trial_outcome,cond_shift}=tmp_stim_col_4(tmp_stim_col_4(:,cond)~=0,cond);
                    Stim_Col_rt_4{trial_trial_outcome,cond_shift}=tmp_stim_col_rt_4(tmp_stim_col_rt_4(:,cond)~=0,cond);
                end
                
            end
        end
        
        for cond=[1:4]
            StimT{cond,1}=vertcat(Stim{cond,1},Stim{cond,4});
            StimT{cond,2}=vertcat(Stim{cond,2},Stim{cond,3});
            StimT{cond,3}=vertcat(Stim{cond,5},Stim{cond,8});
            StimT{cond,4}=vertcat(Stim{cond,6},Stim{cond,7});
            
            Stim_rtT{cond,1}=vertcat(Stim_rt{cond,1},Stim_rt{cond,4});
            Stim_rtT{cond,2}=vertcat(Stim_rt{cond,2},Stim_rt{cond,3});
            Stim_rtT{cond,3}=vertcat(Stim_rt{cond,5},Stim_rt{cond,8});
            Stim_rtT{cond,4}=vertcat(Stim_rt{cond,6},Stim_rt{cond,7});
            
            RespT{cond,1}=vertcat(Resp{cond,1},Resp{cond,4});
            RespT{cond,2}=vertcat(Resp{cond,2},Resp{cond,3});
            RespT{cond,3}=vertcat(Resp{cond,5},Resp{cond,8});
            RespT{cond,4}=vertcat(Resp{cond,6},Resp{cond,7});
            
            Resp_rtT{cond,1}=vertcat(Resp_rt{cond,1},Resp_rt{cond,4});
            Resp_rtT{cond,2}=vertcat(Resp_rt{cond,2},Resp_rt{cond,3});
            Resp_rtT{cond,3}=vertcat(Resp_rt{cond,5},Resp_rt{cond,8});
            Resp_rtT{cond,4}=vertcat(Resp_rt{cond,6},Resp_rt{cond,7});
        end
        
        Stim=StimT;
        Stim_rt=Stim_rtT;
        Resp=RespT;
        Resp_rt=Resp_rtT;
        
        for trial_outcome=1:4
            for cond_blk=1:4
                Stims(trial_outcome,cond_blk,run,subject)=nanmean([Stim{trial_outcome,cond_blk}]);
                Stims_no(trial_outcome,cond_blk,run,subject)=length(Stim{trial_outcome,cond_blk});
                Stims_rt(trial_outcome,cond_blk,run,subject)=nanmean(Stim_rt{trial_outcome,cond_blk});
                Rules(trial_outcome,cond_blk,run,subject)=nanmean(Rule{trial_outcome,cond_blk});
                Rules_no(trial_outcome,cond_blk,run,subject)=length(Rule{trial_outcome,cond_blk});
                Rules_rt(trial_outcome,cond_blk,run,subject)=nanmean(Rule_rt{trial_outcome,cond_blk});
            end
        end
    end
end

clearvars -except runs s Stims Stims_no Stims_rt Rules Rules_no Rules_rt
% save('behavioural_results.mat')
%% Plotting per subject
ccc
% load('behavioural_results.mat')
two_handed_subejcts=1:5;
% one_handed_subejcts=[6:26 28:31];
one_handed_subejcts=[23 12 19 20 26];

Thing=Stims_no/0.4;

% trial_outcome=4; %1=stim error; 2=both error; 3=rule error; 4= correct trials
subplot(1,2,1);
c=0;
for trial_outcome=[1:4]
    c=c+1;
    bars(1)=bar(c,nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,one_handed_subejcts),3),2)),'g','barwidth',0.3);
    hold on;
    errorbar(c,nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,one_handed_subejcts),3),2)),nanstd(nanmean(nanmean(Thing(trial_outcome,:,:,one_handed_subejcts),3),2))./sqrt(length(one_handed_subejcts)),'color','k','capsize',0,'linestyle','none');
    
    bars(2)=bar(c+0.3,nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,two_handed_subejcts),3),2)),'r','barwidth',0.3);
    errorbar(c+0.3,nanmean(nanmean(nanmean(Thing(trial_outcome,:,:,two_handed_subejcts),3),2)),nanstd(nanmean(nanmean(Thing(trial_outcome,:,:,two_handed_subejcts),3),2))./sqrt(length(two_handed_subejcts)),'color','k','capsize',0,'linestyle','none');
    
    if ttest2(squeeze(nanmean(nanmean(Thing(trial_outcome,:,:,one_handed_subejcts),3),2)),squeeze(nanmean(nanmean(Thing(trial_outcome,:,:,two_handed_subejcts),3),2)))
        plot(c+0.15,95,'*b')
    end
    
    data1_temp=squeeze(nanmean(nanmean(Thing(trial_outcome,:,:,one_handed_subejcts),3),2));
    data2=squeeze(nanmean(nanmean(Thing(trial_outcome,:,:,two_handed_subejcts),3),2));
    combinations=nchoosek([1:length(two_handed_subejcts)],5);
    for comb=1:size(combinations,1)
        data1_cmb_mean(comb)=mean(data1_temp(combinations(comb,:)));
    end
    if (nanmean(data2)>0 && sum(data1_cmb_mean<nanmean(data2))>0.99*length(combinations)) || (nanmean(data2)<0 && sum(data1_cmb_mean>nanmean(data2))>0.95*length(combinations))
        plot(c+0.15,20,'*m')
    end
    
end
xticks([[1:4]+0.15])
xticklabels({'Stimulus errors','Both errors','Rule errors','Correct','RT'})
xtickangle(45)
ylabel('Percentage of trials (%)')
legend ([bars(1) bars(2)],{'One-handed','Two-handed'},'location','northwest')
box off; grid on;

subplot(1,2,2);
bars(1)=bar(1,nanmean(nanmean(nanmean(Stims_rt(4,:,:,one_handed_subejcts),3),2)),'g','barwidth',0.3);
hold on;
errorbar(1,nanmean(nanmean(nanmean(Stims_rt(4,:,:,one_handed_subejcts),3),2)),nanstd(nanmean(nanmean(Stims_rt(4,:,:,one_handed_subejcts),3),2))./sqrt(length(one_handed_subejcts)),'color','k','capsize',0,'linestyle','none');
bars(2)=bar(1+0.3,nanmean(nanmean(nanmean(Stims_rt(4,:,:,two_handed_subejcts),3),2)),'r','barwidth',0.3);
errorbar(1+0.3,nanmean(nanmean(nanmean(Stims_rt(4,:,:,two_handed_subejcts),3),2)),nanstd(nanmean(nanmean(Stims_rt(4,:,:,two_handed_subejcts),3),2))./sqrt(length(two_handed_subejcts)),'color','k','capsize',0,'linestyle','none');

if ttest2(squeeze(nanmean(nanmean(Stims_rt(trial_outcome,:,:,one_handed_subejcts),3),2)),squeeze(nanmean(nanmean(Stims_rt(trial_outcome,:,:,two_handed_subejcts),3),2)))
    plot(1+0.15,2.1,'*b')
end

data1_temp=squeeze(nanmean(nanmean(Stims_rt(trial_outcome,:,:,one_handed_subejcts),3),2));
data2=squeeze(nanmean(nanmean(Stims_rt(trial_outcome,:,:,two_handed_subejcts),3),2));
combinations=nchoosek([1:length(two_handed_subejcts)],5);
for comb=1:size(combinations,1)
    data1_cmb_mean(comb)=mean(data1_temp(combinations(comb,:)));
end
if (nanmean(data2)>0 && sum(data1_cmb_mean<nanmean(data2))>0.95*length(combinations)) || (nanmean(data2)<0 && sum(data1_cmb_mean>nanmean(data2))>0.95*length(combinations))
    plot(1+0.15,2.05,'*m')
end

xticks([1+0.15])
xticklabels({'Reaction time'})
xtickangle(45)
ylabel('Reaction time (s)')
legend ([bars(1) bars(2)],{'One-handed','Two-handed'},'location','northwest')
box off; grid on;

