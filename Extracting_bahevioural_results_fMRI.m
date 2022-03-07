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
    tmp_stimT=[];
    tmp_stimT_rt=[];
    s=s+1;
    for run=runs
        clearvars -except StimsAll StimsAll_RT tmp_stimT tmp_stimT_rt runs s Stims Stims_no Stims_rt Rules Rules_no Rules_rt csub Behavioural_directory collapsed rule_version Subjects subject Main_analysis_directory run performances
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

        response_trial_outcome_blk=response_trial_outcome;
        tmp_stim=zeros(size(response_trial_outcome_blk,1),4);
        tmp_stim_rt=nan(size(response_trial_outcome_blk,1),4);
            for trial_trial_outcome=1:4
                % Row1: Stim error;
                % Row2: Stim and Rule error;
                % Row3: Rule error;
                % Row4: Correct;
                Trials_outcome{trial_trial_outcome}=response_trial_outcome_blk(response_trial_outcome_blk(:,trial_trial_outcome)~=0,trial_trial_outcome);
                target_trials=Trials_outcome{trial_trial_outcome}';                               
                
                for i=target_trials
                    tmp_stim(i,trial_trial_outcome)=1;
                    tmp_stim_rt(i,trial_trial_outcome)=test_res.allresults(i,15);                    
                end 
            end
            tmp_stimT=vertcat(tmp_stimT,tmp_stim);
            tmp_stimT_rt=vertcat(tmp_stimT_rt,tmp_stim_rt);
    end    
    StimsAll{subject}=tmp_stimT;
    StimsAll_RT{subject}=tmp_stimT_rt;
end

% s=0;
% for Subject=[1:26 28:31]
%     s=s+1;
%    Stims_subj(s,:)=nanmean(StimsAll{Subject});
%    Stims_subj_rt(s,:)=nanmean(StimsAll_RT{Subject});
% end
% save('behavioural_data_fMRI.mat','Stims_subj','Stims_subj_rt')
