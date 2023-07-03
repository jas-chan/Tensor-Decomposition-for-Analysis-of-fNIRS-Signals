% Excluding Irrelevant Components and Determination of TOI and ROI 
% This code is to apply 2-way ANOVA on the decomposed subjects' signatures after
% applying tensor-decomposition.
% The significant bases will be saved to ANOVA/significant_entity,
% significant_action, significant_interaction
% in this version, ANOVA function was used for unbalanced conditions.
% It will also exclude components on set parameters and combine all the
% components that are left in the according signifcant effects.


%% Data directory and hyperparameters
clc; clear all; % clear command window and clear all variables
for a = 1:10  % goes through all the runs of tensor decomposition
cpd_dir=strcat('./Results/', strcat('CPD_decomp_',string(a))); %the directory of the decompostion after running CPD_decompV4_4.m
load(strcat(cpd_dir,'/CPD_all'), 'U_left','U_right','condition_ind','tensor_right', 'tensor_left');
noSamples=size(U_left{1},1); % number of temporal samples from the spectrogram
no_freq_samp=size(U_left{2},1); % number of spectral samples from the spectrogram
fs=50; % sampling frequency of the fNIRS machine 
window_size=3*fs; % size of the spectrogram window 
overlap=window_size*.9;  % how much to overlap the windows for the spectrgram
signal_length=1351; % number of time points in the original fNIRS data
nfft=2^12; % number of DFT points
time_sp=0:(window_size-overlap)/fs:(signal_length-window_size)/fs; % time (in seconds) of the temporal points from the spectrogram
freq_sp=0:fs/nfft:1;% only the frequencies between 0 and 1

tick_time=(window_size-overlap)/fs; % number of ticks for the temporal mode
tick_freq=nfft/fs; % number of ticks for spectral mode

conditions_title={'Human hand, function', 'Human hand, non-function', 'Mechanical hand, function', ...
    'Mechanical hand, non-function'}; % names of the conditions


%% Finding ANOVA for the signatures 
r_components_l=size(U_left{1},2); % getting the number of components for the left hemisphere
r_components_r=size(U_right{1},2); % getting the number of components for the right hemisphere

% creating tables to store the F-values and p-values 
T_anova_left = table('Size',[r_components_l 6],'VariableNames',{'F_entity', ...
    'P_entity','F_action','P_action','F_interaction','P_interaction'}, ...
    'VariableTypes',{'double', 'double', 'double', 'double', 'double', 'double'});
T_anova_right = table('Size',[r_components_r 6],'VariableNames',{'F_entity', ...
    'P_entity','F_action','P_action','F_interaction','P_interaction'}, ...
    'VariableTypes',{'double', 'double', 'double', 'double', 'double', 'double'});

%make two vectors for each factor of anova
entity_factor=condition_ind;
entity_factor(or(condition_ind==1,condition_ind==2))=1; 
entity_factor(or(condition_ind==3,condition_ind==4))=2; 
action_factor=condition_ind;
action_factor(or(condition_ind==1,condition_ind==3))=1;
action_factor(or(condition_ind==2,condition_ind==4))=2;

%running anova for left hemisphere:
for r=1:r_components_l
    [p,F_table]=anovan(U_left{4}(:,r),{entity_factor,action_factor},...
        'model','interaction','varnames',{'entity','action'},'display','off');
    T_anova_left{r,1}=F_table{2,6}; T_anova_left{r,2}=p(1);
    T_anova_left{r,3}=F_table{3,6}; T_anova_left{r,4}=p(2);
    T_anova_left{r,5}=F_table{4,6}; T_anova_left{r,6}=p(3);
end
%running anova for right hemisphere
for r=1:r_components_r
    [p,F_table]=anovan(U_right{4}(:,r),{entity_factor,action_factor},...
        'model','interaction','varnames',{'entity','action'},'display','off');
    T_anova_right{r,1}=F_table{2,6}; T_anova_right{r,2}=p(1);
    T_anova_right{r,3}=F_table{3,6}; T_anova_right{r,4}=p(2);
    T_anova_right{r,5}=F_table{4,6}; T_anova_right{r,6}=p(3);
end

% Copying the significant patterns or basis to ANOVA/significant_entity,
% significant_action, significant_interaction
entity_dir=strcat(cpd_dir,'/anova','/significant_entity/');
action_dir=strcat(cpd_dir,'/anova','/significant_action/');
inter_dir=strcat(cpd_dir,'/anova','/significant_interaction/');
if ~exist(strcat(cpd_dir,'/anova'))
    mkdir(strcat(cpd_dir,'/anova'))
    mkdir(entity_dir)
    mkdir(action_dir)
    mkdir(inter_dir)
end
img_dir=strcat(cpd_dir,'/figures/');


%% setting exclusion parameters  
noSamples_time=length(U_left{1}); % number of time samples
noSamples_freq=length(U_left{2}); % number of frequency samples
time_exclusion = 15; %time point of end of stimulus presentation (in seconds)
time_exclu_adjust =  find(time_sp==time_exclusion); % changing/adjusting seconds to window sizes to get the right time point
epoch_exclude = .01; % excludes participants if the mean magnitude during stimulus presentation is below this values 
freq_exclude =.1; % the frequency (in Hz) that you want to use for low pass filter (i.e., any components with frequency above this is excluded from analysis)
freq_exclu_adjust = floor(freq_exclude*tick_freq)+1; % ajust the frequency to match the frequency samples
start_time_exclu =0;  % start of stimulus presentation (in seconds)
start_time_exclu_adjust=find(time_sp==start_time_exclu); % changing/adjusting seconds to window sizes to get the right time point

%selecting component that are statistically significant and pass the
%exclusion criteria for the left hemisphere
for r=1:r_components_l
    if T_anova_left{r,2}<0.05 &&(mean(U_left{1}(start_time_exclu_adjust:time_exclu_adjust,r))>epoch_exclude)==1 ...
            && find(ismember(U_left{2}(:,r),max(U_left{2}(:,r))))<=freq_exclu_adjust
        copyfile(strcat(img_dir,'left',num2str(r),'.png'),strcat(entity_dir,'left',num2str(r),'.png'))
    end
    
    if T_anova_left{r,4}<0.05 &&(mean(U_left{1}(start_time_exclu_adjust:time_exclu_adjust,r))>epoch_exclude)==1 ...
            && find(ismember(U_left{2}(:,r),max(U_left{2}(:,r))))<=freq_exclu_adjust
        copyfile(strcat(img_dir,'left',num2str(r),'.png'),strcat(action_dir,'left',num2str(r),'.png'))
    end
    
    if T_anova_left{r,6}<0.05 &&(mean(U_left{1}(start_time_exclu_adjust:time_exclu_adjust,r))>epoch_exclude)==1 ...
            && find(ismember(U_left{2}(:,r),max(U_left{2}(:,r))))<=freq_exclu_adjust
        copyfile(strcat(img_dir,'left',num2str(r),'.png'),strcat(inter_dir,'left',num2str(r),'.png'))
    end
end

%selecting component that are statistically significant and pass the
%exclusion criteria for the left hemisphere
for r=1:r_components_r
    if T_anova_right{r,2}<0.05 &&(mean(U_right{1}(start_time_exclu_adjust:time_exclu_adjust,r))>epoch_exclude)==1 ...
            && find(ismember(U_right{2}(:,r),max(U_right{2}(:,r))))<=freq_exclu_adjust
        copyfile(strcat(img_dir,'right',num2str(r),'.png'),strcat(entity_dir,'right',num2str(r),'.png'))
    end
    
    if T_anova_right{r,4}<0.05 &&(mean(U_right{1}(start_time_exclu_adjust:time_exclu_adjust,r))>epoch_exclude)==1 ...
            && find(ismember(U_right{2}(:,r),max(U_right{2}(:,r))))<=freq_exclu_adjust
        copyfile(strcat(img_dir,'right',num2str(r),'.png'),strcat(action_dir,'right',num2str(r),'.png'))
    end
    
    if T_anova_right{r,6}<0.05 &&(mean(U_right{1}(start_time_exclu_adjust:time_exclu_adjust,r))>epoch_exclude)==1 ...
            && find(ismember(U_right{2}(:,r),max(U_right{2}(:,r))))<=freq_exclu_adjust
        copyfile(strcat(img_dir,'right',num2str(r),'.png'),strcat(inter_dir,'right',num2str(r),'.png'))
    end
end
save(strcat(cpd_dir,'/anova/anova_results'),'T_anova_left','T_anova_right');

end

 
%% summing significant componenent 
cpd_dir='./Results/CPD_decomp_'; %the directory of the decompostion after running tensor_construction_tensor_decomposition.m

%directory to save the final results for each effect
entity_dir=strcat(cpd_dir,'1/anova','/significant_entity/');
action_dir=strcat(cpd_dir,'1/anova','/significant_action/');
inter_dir=strcat(cpd_dir,'1/anova','/significant_interaction/');
save_dir=strcat(cpd_dir,'1/anova');

%main effect of entity
time_entity=zeros(81,1); % time profile for main effect of entity
freq_entity=zeros(82,1); % spectral profile 
channel_entity=zeros(10,1); % ROI 
subj_entity=zeros(70,1); % subject profile 
c_entity=0; % count of the significant main effect components 

%main effect of action
time_action=zeros(81,1);
freq_action=zeros(82,1);
channel_action=zeros(10,1);
subj_action=zeros(70,1);
c_action=0;

%interaction effect
time_inter=zeros(81,1);
freq_inter=zeros(82,1);
channel_inter=zeros(10,1);
subj_inter=zeros(70,1);
c_inter=0;
    
%creates an array of the signifcant components from each run
c_entity_array_l(1,1:10)=0; 
c_action_array_l(1,1:10)=0;
c_inter_array_l(1,1:10)=0;

% left hemisphere
for rep=1:10 %number of runs
    cpd_dir_r=strcat(cpd_dir,num2str(rep)); %directory of the tensor decomposition
    load(strcat(cpd_dir_r,'/CPD_all'), 'U_left','U_right','condition_ind','tensor_right', 'tensor_left'); % loading  tensor decomposition
    load(strcat(cpd_dir_r,'/anova/anova_results'), 'T_anova_left','T_anova_right'); % loading the ANOVA results
    r_components_l=size(U_left{1},2); % loading the number of components

    % keep count of the number of sig components used in one run
    c_entity_run=0;
    c_action_run=0;
    c_inter_run=0;
    

    for r=1:r_components_l % going through all of the components

        %significant main effect of entity type
         if T_anova_left{r,2}<0.05  &&(mean(abs(U_left{1}(1:time_exclu_adjust,r)))>epoch_exclude)==1 ...
            && find(ismember(U_left{2}(:,r),max(U_left{2}(:,r))))<=freq_exclu_adjust %see if component passed the statistically signifcant and exclusion criteria 
            c_entity=c_entity+1; %increase count 
            time_entity(:,c_entity)=U_left{1}(:,r); %load the temporal component
            freq_entity(:,c_entity)=U_left{2}(:,r); %load the spectral component
            channel_entity(:,c_entity)=U_left{3}(:,r); %load the spatial component 
            subj_entity(:,c_entity)=U_left{4}(:,r); %load the subject component
             c_entity_run = c_entity_run+1; %the number of sig components used in one run
        end

        %significant main effect action 
        if T_anova_left{r,4}<0.05 &&(mean(abs(U_left{1}(1:time_exclu_adjust,r)))>epoch_exclude)==1 ...
            && find(ismember(U_left{2}(:,r),max(U_left{2}(:,r))))<=freq_exclu_adjust
           c_action=c_action+1;time_action(:,c_action)=U_left{1}(:,r);
           freq_action(:,c_action)=U_left{2}(:,r);
           channel_action(:,c_action)=U_left{3}(:,r);
           subj_action(:,c_action)=U_left{4}(:,r);
            c_action_run = c_action_run+1;  
        end

        %significant interaction 
        if T_anova_left{r,6}<0.05&&(mean(abs(U_left{1}(1:time_exclu_adjust,r)))>epoch_exclude)==1 ...
            && find(ismember(U_left{2}(:,r),max(U_left{2}(:,r))))<=freq_exclu_adjust
           c_inter=c_inter+1;time_inter(:,c_inter)=U_left{1}(:,r);
           freq_inter(:,c_inter)=U_left{2}(:,r);
           channel_inter(:,c_inter)=U_left{3}(:,r);
           subj_inter(:,c_inter)=U_left{4}(:,r);
           c_inter_run =  c_inter_run +1; 
        end

        %creates an array of the components used in analysis from each run
         c_entity_array_l(1,rep) = c_entity_run; 
         c_action_array_l(1,rep) = c_action_run; 
         c_inter_array_l(1,rep) = c_inter_run;  
    end

end

%mean number of components used for analysis across the ten runs
mean_c_entity_l=mean(c_entity_array_l);
mean_c_action_l=mean(c_action_array_l);
mean_c_inter_l=mean(c_inter_array_l);

%standard deviation of the number of components used for analysis across
%the ten runs 
std_c_entity_l=std(c_entity_array_l);
std_c_action_l=std(c_action_array_l);
std_c_inter_l=std(c_inter_array_l);

%takes the correlation and p values of the components 
[R_time_entity_l, P_time_entity_l] = corrcoef(time_entity);
[R_time_action_l, P_time_action_l] = corrcoef(time_action);
[R_time_inter_l, P_time_inter_l] = corrcoef(time_inter);
[R_channel_entity_l, P_channel_entity_l] = corrcoef(channel_entity);
[R_channel_action_l, P_channel_action_l] = corrcoef(channel_action);
[R_channel_inter_l, P_channel_inter_l] = corrcoef(channel_inter);

%select the temporal and spatial components that are highly correlated (.5)
%and are statistically significant.
selected_temp_entity_l= P_time_entity_l<0.05 * R_time_entity_l>0.5;
selected_temp_action_l= P_time_action_l<0.05 * R_time_action_l>0.5;
selected_temp_inter_l= P_time_inter_l<0.05 * R_time_inter_l>0.5;
selected_channel_entity_l= P_channel_entity_l<0.05 * R_channel_entity_l>0.5;
selected_channel_action_l= P_channel_action_l<0.05 * R_channel_action_l>0.5;
selected_channel_inter_l= P_channel_inter_l<0.05 * R_channel_inter_l>0.5;

% assigning weighting to the highly correlated, statistically significant
% components in proportion of the number of times a similar component has
% occured in other runs (e.g., there are 5 stat sig highly correlated
% components. The three components that are more similar are given more weigting (60%)
% and the other two components that are more similar with each other are
% weighted (40%) for the final results presented. 
select_final_entity_l = selected_temp_entity_l.*selected_channel_entity_l;  
select_final_action_l = selected_temp_action_l.*selected_channel_action_l;
select_final_inter_l = selected_temp_inter_l.*selected_channel_inter_l;
select_final_comp_entity_l = sum(select_final_entity_l, 2)/sum(sum(select_final_entity_l, 2)); 
select_final_comp_action_l = sum(select_final_action_l, 2)/sum(sum(select_final_action_l, 2));
select_final_comp_inter_l = sum(select_final_inter_l, 2)/sum(sum(select_final_inter_l, 2));

%multiplied out the weights for each component (e.g., select_final_comp...) and then summed up the highly correlated components
weighted_time_entity_l = time_entity*select_final_comp_entity_l;
weighted_channel_entity_l = channel_entity*select_final_comp_entity_l;
weighted_freq_entity_l = freq_entity*select_final_comp_entity_l;
weighted_subj_entity_l = subj_entity*select_final_comp_entity_l;

weighted_time_action_l= time_action*select_final_comp_action_l;
weighted_channel_action_l = channel_action*select_final_comp_action_l;
weighted_freq_action_l = freq_action*select_final_comp_action_l;
weighted_subj_action_l = subj_action*select_final_comp_action_l;

weighted_time_inter_l = time_inter*select_final_comp_inter_l;
weighted_channel_inter_l = channel_inter*select_final_comp_inter_l;
weighted_freq_inter_l = freq_inter*select_final_comp_inter_l;
weighted_subj_inter_l = subj_inter*select_final_comp_inter_l;

% plotting the ROI, TOI, spectral profile, and subject profile of the
% significant effects
if weighted_channel_entity_l>0 
    plot_channels(weighted_channel_entity_l,true,entity_dir);
    plot_time(weighted_time_entity_l,true,entity_dir,tick_time);
    plot_freq(weighted_freq_entity_l,true,entity_dir,tick_freq);
    plot_subj(weighted_subj_entity_l,true,entity_dir,condition_ind,conditions_title);
end
if weighted_channel_action_l>0 
    plot_channels(weighted_channel_action_l,true,action_dir);
    plot_time(weighted_time_action_l,true,action_dir,tick_time);
    plot_freq(weighted_freq_action_l,true,action_dir,tick_freq);
    plot_subj(weighted_subj_action_l,true,action_dir,condition_ind,conditions_title);
end
if weighted_channel_inter_l>0
    plot_channels(weighted_channel_inter_l,true,inter_dir);
    plot_time(weighted_time_inter_l,true,inter_dir,tick_time);
    plot_freq(weighted_freq_inter_l,true,inter_dir,tick_freq);
    plot_subj(weighted_subj_inter_l,true,inter_dir,condition_ind,conditions_title);
end

% Combining the significant components like above but for the right hemisphere
channel_entity=zeros(10,1);
c_entity=0;
channel_action=zeros(10,1);
c_action=0;
channel_inter=zeros(10,1);
c_inter=0;

c_entity_array_r(1,1:10)=0; 
c_action_array_r(1,1:10)=0;
c_inter_array_r(1,1:10)=0;

for rep=1:10
    cpd_dir_r=strcat(cpd_dir,num2str(rep));
    load(strcat(cpd_dir_r,'/CPD_all'), 'U_left','U_right','condition_ind','tensor_right', 'tensor_left');
    load(strcat(cpd_dir_r,'/anova/anova_results'), 'T_anova_left','T_anova_right');
    r_components_r=size(U_right{1},2);
       c_entity_run_r=0;
    c_action_run_r=0;
    c_inter_run_r=0;
    for r=1:r_components_r
        if T_anova_right{r,2}<0.05&&(mean(abs(U_right{1}(1:time_exclu_adjust,r)))>epoch_exclude)==1 ...
            && find(ismember(U_right{2}(:,r),max(U_right{2}(:,r))))<=freq_exclu_adjust
            c_entity=c_entity+1; time_entity(:,c_entity)=U_right{1}(:,r);
            freq_entity(:,c_entity)=U_right{2}(:,r);
            channel_entity(:,c_entity)=U_right{3}(:,r);
            subj_entity(:,c_entity)=U_right{4}(:,r);
            c_entity_run_r = c_entity_run_r+1;
        end

        if T_anova_right{r,4}<0.05 &&(mean(abs(U_right{1}(1:time_exclu_adjust,r)))>epoch_exclude)==1 ...
            && find(ismember(U_right{2}(:,r),max(U_right{2}(:,r))))<=freq_exclu_adjust
           c_action=c_action+1; time_action(:,c_action)=U_right{1}(:,r);
           freq_action(:,c_action)=U_right{2}(:,r);
           channel_action(:,c_action)=U_right{3}(:,r);
           subj_action(:,c_action)=U_right{4}(:,r);
           c_action_run_r = c_action_run_r+1;
        end

        if T_anova_right{r,6}<0.05 &&(mean(abs(U_right{1}(1:time_exclu_adjust,r)))>epoch_exclude)==1 ...
            && find(ismember(U_right{2}(:,r),max(U_right{2}(:,r))))<=freq_exclu_adjust
           c_inter=c_inter+1;time_inter(:,c_inter)=U_right{1}(:,r);
           freq_inter(:,c_inter)=U_right{2}(:,r);
           channel_inter(:,c_inter)=U_right{3}(:,r);
           subj_inter(:,c_inter)=U_right{4}(:,r);
           c_inter_run_r =  c_inter_run_r +1; 
        end
        c_entity_array_r(1,rep) = c_entity_run_r;  
        c_action_array_r(1,rep) = c_action_run_r;  
        c_inter_array_r(1,rep) = c_inter_run_r;  
    end
end

%mean number of components used for analysis across the ten runs
mean_c_entity_r=mean(c_entity_array_r);
mean_c_action_r=mean(c_action_array_r);
mean_c_inter_r=mean(c_inter_array_r);

%standard deviation of the number of components used for analysis across
%the ten runs 
std_c_entity_r=std(c_entity_array_r);
std_c_action_r=std(c_action_array_r);
std_c_inter_r=std(c_inter_array_r);

%takes the correlation and p values of the components 
[R_time_entity_r, P_time_entity_r] = corrcoef(time_entity);
[R_time_action_r, P_time_action_r] = corrcoef(time_action);
[R_time_inter_r, P_time_inter_r] = corrcoef(time_inter);
[R_channel_entity_r, P_channel_entity_r] = corrcoef(channel_entity);
[R_channel_action_r, P_channel_action_r] = corrcoef(channel_action);
[R_channel_inter_r, P_channel_inter_r] = corrcoef(channel_inter);

%select the temporal and spatial components that are highly correlated (.5)
%and are statistically significant.
selected_temp_entity_r= P_time_entity_r<0.05 * R_time_entity_r>0.5;
selected_temp_action_r= P_time_action_r<0.05 * R_time_action_r>0.5;
selected_temp_inter_r= P_time_inter_r<0.05 * R_time_inter_r>0.5;
selected_channel_entity_r= P_channel_entity_r<0.05 * R_channel_entity_r>0.5;
selected_channel_action_r= P_channel_action_r<0.05 * R_channel_action_r>0.5;
selected_channel_inter_r= P_channel_inter_r<0.05 * R_channel_inter_r>0.5;
select_final_entity_r = selected_temp_entity_r.*selected_channel_entity_r;
select_final_action_r = selected_temp_action_r.*selected_channel_action_r;
select_final_inter_r = selected_temp_inter_r.*selected_channel_inter_r;
select_final_comp_entity_r = sum(select_final_entity_r, 2)/sum(sum(select_final_entity_r, 2)); 
select_final_comp_action_r = sum(select_final_action_r, 2)/sum(sum(select_final_action_r, 2));
select_final_comp_inter_r = sum(select_final_inter_r, 2)/sum(sum(select_final_inter_r, 2));

%multiplied out the weights for each component (e.g., select_final_comp...) and then summed up the highly correlated components
weighted_time_entity_r = time_entity*select_final_comp_entity_r;
weighted_channel_entity_r = channel_entity*select_final_comp_entity_r;
weighted_freq_entity_r = freq_entity*select_final_comp_entity_r;
weighted_subj_entity_r = subj_entity*select_final_comp_entity_r;
weighted_time_action_r= time_action*select_final_comp_action_r;
weighted_channel_action_r = channel_action*select_final_comp_action_r;
weighted_freq_action_r = freq_action*select_final_comp_action_r;
weighted_subj_action_r = subj_action*select_final_comp_action_r;
weighted_time_inter_r = time_inter*select_final_comp_inter_r;
weighted_channel_inter_r = channel_inter*select_final_comp_inter_r;
weighted_freq_inter_r = freq_inter*select_final_comp_inter_r;
weighted_subj_inter_r = subj_inter*select_final_comp_inter_r;

if weighted_channel_entity_r>0 
    plot_channels(weighted_channel_entity_r,false,entity_dir);
    plot_time(weighted_time_entity_r,false,entity_dir,tick_time);
    plot_freq(weighted_freq_entity_r,false,entity_dir,tick_freq);
    plot_subj(weighted_subj_entity_r,false,entity_dir,condition_ind,conditions_title);
    
end
if weighted_channel_action_r>0 
    plot_channels(weighted_channel_action_r,false,action_dir);
    plot_time(weighted_time_action_r,false,action_dir,tick_time);
    plot_freq(weighted_freq_action_r,false,action_dir,tick_freq);
    plot_subj(weighted_subj_action_r,false,action_dir,condition_ind,conditions_title);
end
if weighted_channel_inter_r>0
    plot_channels(weighted_channel_inter_r,false,inter_dir);
    plot_time(weighted_time_inter_r,false,inter_dir,tick_time);
    plot_freq(weighted_freq_inter_r,false,inter_dir,tick_freq);
    plot_subj(weighted_subj_inter_r,false,inter_dir,condition_ind,conditions_title);
end


save(fullfile(save_dir,'components_analysis'),'mean_c_entity_l',...
    'mean_c_action_l','mean_c_inter_l','std_c_entity_l','std_c_action_l',...
    'std_c_inter_l','mean_c_entity_r','mean_c_action_r','mean_c_inter_r',...
    'std_c_entity_r','std_c_action_r','std_c_inter_r',...
    'select_final_comp_entity_l','select_final_comp_action_l',...
    'select_final_comp_inter_l','select_final_comp_entity_r',... 
    'select_final_comp_action_r','select_final_comp_inter_r',...
    'weighted_time_entity_l','weighted_time_action_l','weighted_time_inter_l',...
    'weighted_time_entity_r','weighted_time_action_r','weighted_time_inter_r',...
    'noSamples_time','tick_time'); % saving variables 

% plotting ROI
function plot_channels(data,left,saveDir)
    FigH = figure('Position', get(0, 'Screensize'),'visible',false);
    xbase=1;ybase=1; activations=zeros(4+xbase,8+ybase);
    Infimg_left=imread('.\left pic.gif',1);
    Infimg_left=imresize(Infimg_left,0.2); 
    Infimg_right=imread('.\right pic.gif',1);
    Infimg_right=imresize(Infimg_right,0.5); 
    if left==true
        activations(xbase+3,ybase+[2,4,6])=data(1:3); 
        activations(xbase+2,ybase+[1,3,5,7])=data(4:7);
        activations(xbase+1,ybase+[2,4,6])=data(8:10);

        %ploting the emitters and detectors placed on left hemisphere
        plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ks','MarkerFaceColor','k','MarkerSize', 10)
        hold on
        plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ro','MarkerFaceColor','r','MarkerSize', 10)

        x = [1:8+xbase];
        y = [1:4+ybase];  

        levels=linspace(min(min(activations)),max(max(activations)),10); 
        contour(x,y,activations,levels,'ShowText','off','LineWidth',7) 
        contourcmap('hsv', levels,'Colorbar', 'on', ...
        'Location', 'horizontal', ...
        'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(data)))...
            ,' to ',num2str(max(max(data)))));
        axis equal                                
        axis square off
        imgh=imshow(repmat(Infimg_left,1,1,3), 'XData', [-11 14], 'YData', [-9 15]);
        imgh.AlphaData = .4;
        hold off
        xlim([0.77 9.716])
        ylim([0.8543 7.795])
    else
        activations(xbase+1,ybase+[6,4,2])=data(8:10);
        activations(xbase+2,ybase+[7,5,3,1])=data(4:7);
        activations(xbase+3,ybase+[6,4,2])=data(1:3);  

        %ploting the emitters and detectors placed on left hemisphere
        plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ro','MarkerFaceColor','k','MarkerSize', 10)
        hold on
        plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ks','MarkerFaceColor','r','MarkerSize', 10)

        x = [1:8+xbase];
        y = [1:4+ybase];  

        levels=linspace(min(min(activations)),max(max(activations)),10);  
        contour(x,y,activations,levels,'ShowText','off','LineWidth',7)  
        contourcmap('hsv', levels,'Colorbar', 'on', ...
       'Location', 'horizontal'); 
        axis equal                                 
        axis square off
        imgh=imshow(repmat(Infimg_right,1,1,3), 'XData', [-4 21], 'YData', [-9 15]);
        imgh.AlphaData = .4;
        hold off
        xlim([0.77 9.716])
        ylim([0.8543 7.795]) 
    end
    
    F    = getframe(FigH);
    if left
        imwrite(F.cdata, fullfile(saveDir,'channels_l_trials_selected_components.png'), 'png')
        saveas(FigH,fullfile(saveDir,'channels_l_trials_selected_components.fig'), 'fig')
        close 
    else
        imwrite(F.cdata, fullfile(saveDir,'channels_r_trials_selected_components.png'), 'png')
        saveas(FigH,fullfile(saveDir,'channels_r_trials_selected_components.fig'), 'fig')
        close 
    end
end


%plotting temporal profile
function plot_time(data,left,saveDir,tick_time)
    noSamples_time=length(data);
    FigH = figure('Position', get(0, 'Screensize'),'visible',false);
    plot(data,'LineWidth',8,'color','black')

    xticks(0:1/tick_time:noSamples_time)
    xticklabels(0:noSamples_time*tick_time)
    set(gca,'box','off', 'FontSize', 30);
    F    = getframe(FigH);
    if left
        imwrite(F.cdata, fullfile(saveDir,'time_l_trials_selected_components.png'), 'png')
        saveas(FigH,fullfile(saveDir,'time_l_trials_selected_components.fig'), 'fig')
        close 
    else
        imwrite(F.cdata, fullfile(saveDir,'time_r_trials_selected_components.png'), 'png')
         saveas(FigH,fullfile(saveDir,'time_r_trials_selected_components.fig'), 'fig')
        close 
    end
end

% plotting spectral profile
function plot_freq(data,left,saveDir,tick_freq)
    noSamples_freq=length(data);
    FigH = figure('Position', get(0, 'Screensize'),'visible',false);
    plot(data,'LineWidth',8,'color','black')
    xticks(0:tick_freq:noSamples_freq)
    xticklabels(0:noSamples_freq/tick_freq)
    set(gca,'box','off', 'FontSize', 30);
    F    = getframe(FigH);
    if left
        imwrite(F.cdata, fullfile(saveDir,'freq_l_trials_selected_components.png'), 'png')
        saveas(FigH,fullfile(saveDir,'freq_l_trials_selected_components.fig'), 'fig')
        close 
    else
        imwrite(F.cdata, fullfile(saveDir,'freq_r_trials_selected_components.png'), 'png')
        saveas(FigH,fullfile(saveDir,'freq_r_trials_selected_components.fig'), 'fig')
        close 
    end
end

% plotting subject profile 
function plot_subj(data,left,saveDir,condition_ind,conditions) 
    FigH = figure('Position', get(0, 'Screensize'),'visible',false);
    plot(data,'LineWidth',8)
    starti=1;
    for conditionI=1:length(conditions)
        endi=starti+length(find(condition_ind==conditionI))-1;
        bar(conditionI,mean(data(starti:endi)));
        hold on;
        starti=endi+1;
    end    

     xlabel('Condition','Fontsize',18,'FontName','Times New Roman')
     set(gca,'box','off', 'FontSize', 30, 'XTickLabel',[]);
     legend(conditions,'Location','best');
        
     F    = getframe(FigH);
    if left
         imwrite(F.cdata, fullfile(saveDir,'subj_l_trials_selected_components.png'), 'png')
         saveas(FigH,fullfile(saveDir,'subj_l_trials_selected_components.fig'), 'fig')
         close 
    else
         imwrite(F.cdata, fullfile(saveDir,'subj_r_trials_selected_components.png'), 'png')
         saveas(FigH, fullfile(saveDir,'subj_r_trials_selected_components.fig'), 'fig')
         close 
    end 
end



