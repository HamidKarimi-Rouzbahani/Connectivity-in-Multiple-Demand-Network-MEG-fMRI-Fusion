clc;
close all;
clear all;
Main_analysis_directory='D:/Hamid/Postdoc/MD/Analyses';
subject=1;
Behavioural_directory='D:\Hamid\Postdoc\MD\fMRI_Behavioural_data\';


runss=dir([Behavioural_directory,'/test_res_p',num2str(subject),'_r*.mat']);
runs=[1:size(runss,1)];

rule_version=[1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2];

% collapsed=1;

for run=runs
    clearvars -except csub Behavioural_directory collapsed rule_version Subjects subject Main_analysis_directory run performances
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
    response_type=zeros(length(test_res.allresults(:,12)),4);
    % Trial type categorization
    tmp_ro=zeros(length(test_res.allresults(:,12)),4);
    
    for i=all_responded_trials
        
        if test_res.allresults(i,4)==rule_version(subject)
            
            if test_res.allresults(i,7)==1 && test_res.allresults(i,9)==1
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==2
                % Stim + Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==3
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
            elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==4
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
                
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==1
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==2
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==3
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==4
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
                
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==1
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==2
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==3
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==4
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
                
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==1
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==2
                % Corrrect
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==3
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==4
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            end
            
        else
            if test_res.allresults(i,7)==1 && test_res.allresults(i,9)==1
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==2
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==3
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
            elseif test_res.allresults(i,7)==1 && test_res.allresults(i,9)==4
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
                
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==1
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==2
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==3
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            elseif test_res.allresults(i,7)==2 && test_res.allresults(i,9)==4
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
                
                
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==1
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==2
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==3
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
            elseif test_res.allresults(i,7)==3 && test_res.allresults(i,9)==4
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
                
                
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==1
                % Correct
                tmp_ro(i,4)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,4)=i;
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==2
                % Rule
                tmp_ro(i,3)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,3)=i;
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==3
                % Stim+Rule
                tmp_ro(i,2)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,2)=i;
            elseif test_res.allresults(i,7)==4 && test_res.allresults(i,9)==4
                % Stim
                tmp_ro(i,1)=test_res.allresults(i,12)+first_stim_delay;
                response_type(i,1)=i;
            end
        end
        %         Responded_RT_tmp(i,find(tmp_ro(i,:)~=0))=test_res.allresults(i,15);
    end
    
    for blk=1:2
        response_type_blk=response_type([[1:80]+(blk-1)*80],:);
        
        for trial_type=1:4
            % Row1: Stim error;
            % Row2: Stim and Rule error;
            % Row3: Rule error;
            % Row4: Correct;
            Trials_outcome{trial_type}=response_type_blk(response_type_blk(:,trial_type)~=0,trial_type);
            target_trials=Trials_outcome{trial_type}';
            
            
            tmp_stim=zeros(size(response_type_blk,1),4);
            tmp_rule=zeros(size(response_type_blk,1),2);
            tmp_rule_color=zeros(size(response_type_blk,1),4);
            tmp_resp=zeros(size(response_type_blk,1),4);
            tmp_stim_rt=zeros(size(response_type_blk,1),4);
            tmp_rule_rt=zeros(size(response_type_blk,1),2);
            tmp_rule_color_rt=zeros(size(response_type_blk,1),4);
            tmp_resp_rt=zeros(size(response_type_blk,1),4);
            tmp_stim_col_16=zeros(size(response_type_blk,1),16);
            tmp_stim_col_rt_16=zeros(size(response_type_blk,1),16);
            tmp_stim_col_8=zeros(size(response_type_blk,1),8);
            tmp_stim_col_rt_8=zeros(size(response_type_blk,1),8);
            tmp_stim_col_4=zeros(size(response_type_blk,1),4);
            tmp_stim_col_rt_4=zeros(size(response_type_blk,1),4);
            tmp_stim_side_col_4=zeros(size(response_type_blk,1),4);
            tmp_stim_side_col_rt_4=zeros(size(response_type_blk,1),4);
            
            tmp_resp_col_4=zeros(size(response_type_blk,1),4);
            tmp_resp_col_rt_4=zeros(size(response_type_blk,1),4);

            tmp_stim_col_clpsd_4=zeros(size(response_type_blk,1),4);
            tmp_stim_col_clpsd_rt_4=zeros(size(response_type_blk,1),4);
            
            tmp_resp_col_clpsd_4=zeros(size(response_type_blk,1),4);
            tmp_resp_col_clpsd_rt_4=zeros(size(response_type_blk,1),4);
            
            tmp_respInOut_col_clpsd_4=zeros(size(response_type_blk,1),4);
            tmp_respInOut_col_clpsd_rt_4=zeros(size(response_type_blk,1),4);
            
            tmp_respInOut_col_clpsd2_4=zeros(size(response_type_blk,1),4);
            tmp_respInOut_col_clpsd2_rt_4=zeros(size(response_type_blk,1),4);
            
            tmp_stim_side_col_clpsd_4=zeros(size(response_type_blk,1),4);
            tmp_stim_side_col_clpsd_rt_4=zeros(size(response_type_blk,1),4);
            
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

                % corrected
                bias=(test_res.allresults(i,4)-1)*2;
                if test_res.allresults(i,7)==3
                    tmp=2;
                elseif test_res.allresults(i,7)==4
                    tmp=1;
                elseif test_res.allresults(i,7)<3
                    tmp=test_res.allresults(i,7);
                end
                tmp_stim_col_4(t,tmp+bias)=test_res.allresults(i,12)+first_stim_delay;
                tmp_stim_col_rt_4(t,tmp+bias)=test_res.allresults(i,15);

                bias=(test_res.allresults(i,4)-1)*2;
                if test_res.allresults(i,7)>2
                    tmp=2;
                elseif test_res.allresults(i,7)<3
                    tmp=1;
                end
                tmp_stim_side_col_4(t,tmp+bias)=test_res.allresults(i,12)+first_stim_delay;
                tmp_stim_side_col_rt_4(t,tmp+bias)=test_res.allresults(i,15);

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
                
                bias=(test_res.allresults(i,4)-1)*2;
                if test_res.allresults(i,9)==1 || test_res.allresults(i,9)==2
                    tmp=1;
                elseif test_res.allresults(i,9)==3 || test_res.allresults(i,9)==4
                    tmp=2;                
                end
                tmp_resp_col_4(t,tmp+bias)=test_res.allresults(i,12)+first_stim_delay;
                tmp_resp_col_rt_4(t,tmp+bias)=test_res.allresults(i,15);

                
                bias1=(test_res.allresults(i,5)-1)*2; % across rule colors (combining rules)
                if test_res.allresults(i,7)==3
                    tmp=2;
                elseif test_res.allresults(i,7)==4
                    tmp=1;
                elseif test_res.allresults(i,7)<3
                    tmp=test_res.allresults(i,7);
                end
                tmp_stim_col_clpsd_4(t,tmp+bias1)=test_res.allresults(i,12)+first_stim_delay;
                tmp_stim_col_clpsd_rt_4(t,tmp+bias1)=test_res.allresults(i,15);
                           
                if test_res.allresults(i,4)>1 && test_res.allresults(i,5)==2
                    tmp1=1;
                elseif test_res.allresults(i,4)>1 && test_res.allresults(i,5)==1
                    tmp1=2;
                else
                    tmp1=test_res.allresults(i,5);
                end
                bias2=(tmp1-1)*2;
                if test_res.allresults(i,9)==1 || test_res.allresults(i,9)==2
                    tmp2=1;
                elseif test_res.allresults(i,9)==3 || test_res.allresults(i,9)==4
                    tmp2=2;
                end
                tmp_resp_col_clpsd_4(t,tmp2+bias2)=test_res.allresults(i,12)+first_stim_delay;
                tmp_resp_col_clpsd_rt_4(t,tmp2+bias2)=test_res.allresults(i,15);
                
                bias2=(tmp1-1)*2;
                if test_res.allresults(i,9)==1 || test_res.allresults(i,9)==4
                    tmp2=1;
                elseif test_res.allresults(i,9)==3 || test_res.allresults(i,9)==2
                    tmp2=2;
                end
                tmp_respInOut_col_clpsd_4(t,tmp2+bias2)=test_res.allresults(i,12)+first_stim_delay;
                tmp_respInOut_col_clpsd_rt_4(t,tmp2+bias2)=test_res.allresults(i,15);

                
                if test_res.allresults(i,4)>1 && test_res.allresults(i,5)==1
                    tmp1=1;
                elseif test_res.allresults(i,4)>1 && test_res.allresults(i,5)==2
                    tmp1=2;
                else
                    tmp1=test_res.allresults(i,5);
                end
                bias2=(tmp1-1)*2;
                if test_res.allresults(i,9)==1 || test_res.allresults(i,9)==4
                    tmp2=1;
                elseif test_res.allresults(i,9)==3 || test_res.allresults(i,9)==2
                    tmp2=2;
                end
                tmp_respInOut_col_clpsd2_4(t,tmp2+bias2)=test_res.allresults(i,12)+first_stim_delay;
                tmp_respInOut_col_clpsd2_rt_4(t,tmp2+bias2)=test_res.allresults(i,15);
                
                if test_res.allresults(i,7)==1 || test_res.allresults(i,7)==2
                    tmp=1;
                elseif test_res.allresults(i,7)==3 || test_res.allresults(i,7)==4
                    tmp=2;
                end
                tmp_stim_side_col_clpsd_4(t,tmp+bias)=test_res.allresults(i,12)+first_stim_delay;
                tmp_stim_side_col_clpsd_rt_4(t,tmp+bias)=test_res.allresults(i,15);
            end
            for cond=[1:4]
                cond_shift=cond+(blk-1)*4;
                Stim_4_stims_2_Blk{trial_type,cond_shift}=tmp_stim(tmp_stim(:,cond)~=0,cond);
                Stim_rt_4_stims_2_Blk{trial_type,cond_shift}=tmp_stim_rt(tmp_stim_rt(:,cond)~=0,cond);
                
                if cond<3
                    cond_shiftR=cond+(blk-1)*2;
                    Rule{trial_type,cond_shiftR}=tmp_rule(tmp_rule(:,cond)~=0,cond);
                    Rule_rt{trial_type,cond_shiftR}=tmp_rule_rt(tmp_rule_rt(:,cond)~=0,cond);
                end
                Rule_color{trial_type,cond_shift}=tmp_rule_color(tmp_rule_color(:,cond)~=0,cond);
                Rule_color_rt{trial_type,cond_shift}=tmp_rule_color_rt(tmp_rule_color_rt(:,cond)~=0,cond);
                
                Resp{trial_type,cond_shift}=tmp_resp(tmp_resp(:,cond)~=0,cond);
                Resp_rt{trial_type,cond_shift}=tmp_resp_rt(tmp_resp_rt(:,cond)~=0,cond);
            end
            
            for cond=[1:16]
                cond_shift=cond+(blk-1)*16;
                Stim_Col_16{trial_type,cond_shift}=tmp_stim_col_16(tmp_stim_col_16(:,cond)~=0,cond);
                Stim_Col_rt_16{trial_type,cond_shift}=tmp_stim_col_rt_16(tmp_stim_col_rt_16(:,cond)~=0,cond);
            end
            
            for cond=[1:8]
                cond_shift=cond+(blk-1)*8;
                Stim_Col_8{trial_type,cond_shift}=tmp_stim_col_8(tmp_stim_col_8(:,cond)~=0,cond);
                Stim_Col_rt_8{trial_type,cond_shift}=tmp_stim_col_rt_8(tmp_stim_col_rt_8(:,cond)~=0,cond);
            end
            
            for cond=[1:4]
                cond_shift=cond+(blk-1)*4;
                Stim_Col_4{trial_type,cond_shift}=tmp_stim_col_4(tmp_stim_col_4(:,cond)~=0,cond);
                Stim_Col_rt_4{trial_type,cond_shift}=tmp_stim_col_rt_4(tmp_stim_col_rt_4(:,cond)~=0,cond);
                
                Stim_side_Col_4{trial_type,cond_shift}=tmp_stim_side_col_4(tmp_stim_col_4(:,cond)~=0,cond);
                Stim_side_Col_rt_4{trial_type,cond_shift}=tmp_stim_side_col_rt_4(tmp_stim_col_rt_4(:,cond)~=0,cond);
                
                
                Resp_Col_4{trial_type,cond_shift}=tmp_resp_col_4(tmp_resp_col_4(:,cond)~=0,cond);
                Resp_Col_rt_4{trial_type,cond_shift}=tmp_resp_col_rt_4(tmp_resp_col_rt_4(:,cond)~=0,cond);
                
                Stim_Col_Cues_Colpsd_4{trial_type,cond_shift}=tmp_stim_col_clpsd_4(tmp_stim_col_clpsd_4(:,cond)~=0,cond);
                Stim_Col_Cues_Colpsd_rt_4{trial_type,cond_shift}=tmp_stim_col_clpsd_rt_4(tmp_stim_col_clpsd_rt_4(:,cond)~=0,cond);

                Resp_Col_Cues_Colpsd_4{trial_type,cond_shift}=tmp_resp_col_clpsd_4(tmp_resp_col_clpsd_4(:,cond)~=0,cond);
                Resp_Col_Cues_Colpsd_rt_4{trial_type,cond_shift}=tmp_resp_col_clpsd_rt_4(tmp_resp_col_clpsd_rt_4(:,cond)~=0,cond);

              
                Resp_InOut_Cues_Colpsd_4{trial_type,cond_shift}=tmp_respInOut_col_clpsd_4(tmp_respInOut_col_clpsd_4(:,cond)~=0,cond);
                Resp_InOut_Cues_Colpsd_rt_4{trial_type,cond_shift}=tmp_respInOut_col_clpsd_rt_4(tmp_respInOut_col_clpsd_rt_4(:,cond)~=0,cond);
                
                Resp_InOut_Cues_Colpsd2_4{trial_type,cond_shift}=tmp_respInOut_col_clpsd2_4(tmp_respInOut_col_clpsd2_4(:,cond)~=0,cond);
                Resp_InOut_Cues_Colpsd2_rt_4{trial_type,cond_shift}=tmp_respInOut_col_clpsd2_rt_4(tmp_respInOut_col_clpsd2_rt_4(:,cond)~=0,cond);
             
                
                Stim_side_Col_Cues_Colpsd_4{trial_type,cond_shift}=tmp_stim_side_col_clpsd_4(tmp_stim_side_col_clpsd_4(:,cond)~=0,cond);
                Stim_side_Col_Cues_Colpsd_rt_4{trial_type,cond_shift}=tmp_stim_side_col_clpsd_rt_4(tmp_stim_side_col_clpsd_rt_4(:,cond)~=0,cond);

            end
        end
    end
    
    for cond=[1:4]
        StimT{cond,1}=vertcat(Stim_4_stims_2_Blk{cond,1},Stim_4_stims_2_Blk{cond,4});
        StimT{cond,2}=vertcat(Stim_4_stims_2_Blk{cond,2},Stim_4_stims_2_Blk{cond,3});
        StimT{cond,3}=vertcat(Stim_4_stims_2_Blk{cond,5},Stim_4_stims_2_Blk{cond,8});
        StimT{cond,4}=vertcat(Stim_4_stims_2_Blk{cond,6},Stim_4_stims_2_Blk{cond,7});
        
        Stim_rtT{cond,1}=vertcat(Stim_rt_4_stims_2_Blk{cond,1},Stim_rt_4_stims_2_Blk{cond,4});
        Stim_rtT{cond,2}=vertcat(Stim_rt_4_stims_2_Blk{cond,2},Stim_rt_4_stims_2_Blk{cond,3});
        Stim_rtT{cond,3}=vertcat(Stim_rt_4_stims_2_Blk{cond,5},Stim_rt_4_stims_2_Blk{cond,8});
        Stim_rtT{cond,4}=vertcat(Stim_rt_4_stims_2_Blk{cond,6},Stim_rt_4_stims_2_Blk{cond,7});

        StimTside{cond,1}=vertcat(Stim_4_stims_2_Blk{cond,1},Stim_4_stims_2_Blk{cond,2});
        StimTside{cond,2}=vertcat(Stim_4_stims_2_Blk{cond,3},Stim_4_stims_2_Blk{cond,4});
        StimTside{cond,3}=vertcat(Stim_4_stims_2_Blk{cond,5},Stim_4_stims_2_Blk{cond,6});
        StimTside{cond,4}=vertcat(Stim_4_stims_2_Blk{cond,7},Stim_4_stims_2_Blk{cond,8});
        
        Stim_rtTside{cond,1}=vertcat(Stim_rt_4_stims_2_Blk{cond,1},Stim_rt_4_stims_2_Blk{cond,2});
        Stim_rtTside{cond,2}=vertcat(Stim_rt_4_stims_2_Blk{cond,3},Stim_rt_4_stims_2_Blk{cond,4});
        Stim_rtTside{cond,3}=vertcat(Stim_rt_4_stims_2_Blk{cond,5},Stim_rt_4_stims_2_Blk{cond,6});
        Stim_rtTside{cond,4}=vertcat(Stim_rt_4_stims_2_Blk{cond,7},Stim_rt_4_stims_2_Blk{cond,8});
        
        RespTdec{cond,1}=vertcat(Resp{cond,1},Resp{cond,4});
        RespTdec{cond,2}=vertcat(Resp{cond,2},Resp{cond,3});
        RespTdec{cond,3}=vertcat(Resp{cond,5},Resp{cond,8});
        RespTdec{cond,4}=vertcat(Resp{cond,6},Resp{cond,7});
        
        Resp_rtTdec{cond,1}=vertcat(Resp_rt{cond,1},Resp_rt{cond,4});
        Resp_rtTdec{cond,2}=vertcat(Resp_rt{cond,2},Resp_rt{cond,3});
        Resp_rtTdec{cond,3}=vertcat(Resp_rt{cond,5},Resp_rt{cond,8});
        Resp_rtTdec{cond,4}=vertcat(Resp_rt{cond,6},Resp_rt{cond,7});
    end
    Stim=StimT;
    Stim_rt=Stim_rtT;
    Respdec=RespTdec;
    Resp_rtdec=Resp_rtTdec;

    save(fullfile([Main_analysis_directory,'/Results_temp_playing/'],...
        num2str(csub,'S%03d'),['Behav_Subj_',num2str(subject),'_run_',...
        num2str(run),'_per_blk.mat']),'Resp_InOut_Cues_Colpsd2_4','Resp_InOut_Cues_Colpsd2_rt_4','Respdec','Resp_rtdec','StimTside','Stim_rtTside','Resp_InOut_Cues_Colpsd_4','Resp_InOut_Cues_Colpsd_rt_4','Stim_side_Col_4','Stim_side_Col_rt_4','Resp_Col_Cues_Colpsd_4','Resp_Col_Cues_Colpsd_rt_4','Resp_Col_4','Resp_Col_rt_4','Stim_side_Col_Cues_Colpsd_4','Stim_side_Col_Cues_Colpsd_rt_4','Stim_Col_Cues_Colpsd_4','Stim_Col_Cues_Colpsd_rt_4','Stim_4_stims_2_Blk','Stim_rt_4_stims_2_Blk',...
        'Stim','Rule','Rule_color','Resp','Stim_rt','Rule_rt','Rule_color_rt','Resp_rt',...
        'Stim_Col_16','Stim_Col_rt_16','Stim_Col_8','Stim_Col_rt_8','Stim_Col_4','Stim_Col_rt_4');
end

done=1;
