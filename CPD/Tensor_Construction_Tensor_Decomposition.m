% This script is to apply between-subject analysis for all the conditions and 
% for each hemisphere separatly. Using the data filtered using 0.1Hz LPF.
% In this version, time-frequency matrix was added to the tensor.
% This runs a non-negative CPD.
% This has different components for the left and right hemisphere tensors.


%% Data directory and hyperparameters
clc; clear all; % clear command window and clear all variables
data_dir_new='.\Text Files\'; % file path of the fNIRS data
conditions={'Hand function_text','Hand motion_text','Robot function_text','Robot motion_text'}; % names of the folders for each condition
conditions_title={'Human hand, function', 'Human hand, non-function', 'Mechanical claw, function', ...
    'Mechanical claw, non-function'}; % actual names of the conditions
fs=50; % sampling rate of the fNIRS machine
window_size=3*fs; % size of the spectrogram window 
overlap=window_size*.9; % how much to overlap the windows for the spectrgram
signal_length=1351; % number of time points in the original fNIRS data
nfft=2^12; % number of DFT points

r_components_l=50; % number of components for the left hemisphere
r_components_r=50; % number of components for the right hemipshere

for a = 1:10 % number of runs 
saveDir=strcat(fullfile('Results',strcat('CPD_decomp_',string(a)))); % where the decomposition is saved 
if ~exist(saveDir)
    mkdir(saveDir)
end

%% Tensor Construction
tensor=zeros(81,82,20,0); % dimensions of the tensor (temporal by spectral by spatial by subject)
condition_ind=zeros(0); % keeping track of condition number
count_subj=1; % keep track of subject number
for conditionI=1:4 % go through all of the conditions
    tensor_temp=read_condition_tensor(data_dir_new,conditions{conditionI}); % read the fNIRS files 
    tensor_temp = fillmissing(tensor_temp,'constant',0);     %dealing with the missing data (NaN)
    for subjI=1:size(tensor_temp,3) % go through all of the subjects
        for channI=1:20 % goes through all channels
            spec_chann=spectrogram(tensor_temp(:,channI,subjI),kaiser(window_size,5),round(overlap),nfft,fs); % running the spectrogram on the (temporal by spatial by subject) fNIRS data 
            spec_chann =abs(spec_chann(1:ceil(nfft*1/fs),:)); % Using only the frequencies between 0 and 1
            tensor(:,:,channI,count_subj)=spec_chann'; % creates a tensor with temporal by spectral by spatial by subject dimensions
        end
        count_subj=count_subj+1; % increase subject number  
    end
    condition_ind(end+1:end+size(tensor_temp,3))=conditionI; % increase the condition number 
end


%Splitting the tensor into 2, one for each hemisphere
tensor_left=tensor(:,:,1:10,:);
tensor_right=tensor(:,:,11:end,:);

%% Tensor Decomposition 

% %Uncomment the code below to identify the stop criteria to find number of
% %component (r_components) to extract for the left hemisphere by finding the  
% % L-curve corner using the rankest function. 
% options.MaxR=75; % setting the max number to components to extract
% [r_components_l,L_lb_l,L_cpd_l]=rankest(tensor_left,options); %the L-curve ranks L_lb(:,1)
%     %and L_cpd(:,1) and the corresponding lower bound on the truncation
%     %error L_lb(:,2) and relative error of the CPD approximation L_cpd(:,2),
%     %respectively.
% plot(L_cpd_l(:,1),L_cpd_l(:,2)); hold on; plot(L_lb_l(1:options.MaxR,1),L_lb_l(1:options.MaxR,2));
% ylabel('The relative error of the CPD-left hemisphere'); xlabel('The number of rank-one terms R');
% legend('Relative error of the CPD approximation', 'Lower bound on the truncation error')

% %Uncomment the code below to identify the stop criteria to find number of
% %component (r_components) to extract for the right hemisphere by finding the  
% % L-curve corner using the rankest function. 
% options.MaxR=75;
% [r_components_r,L_lb_r,L_cpd_r]=rankest(tensor_left,options); %the L-curve ranks L_lb(:,1)
%     %and L_cpd(:,1) and the corresponding lower bound on the truncation
%     %error L_lb(:,2) and relative error of the CPD approximation L_cpd(:,2),
%     %respectively.
% plot(L_cpd_r(:,1),L_cpd_r(:,2)); hold on; plot(L_lb_r(1:options.MaxR,1),L_lb_r(1:options.MaxR,2));
% ylabel('The relative error of the CPD-left hemisphere'); xlabel('The number of rank-one terms R');
% legend('Relative error of the CPD approximation', 'Lower bound on the truncation error')
% %Stop criteria to find r_components is L-curve corner using rankest


% left hemisphere
Sol_left=non_negative_cpd(tensor_left,r_components_l); % runnning a non-negative CPD
U_left={Sol_left.factors.A,Sol_left.factors.B, Sol_left.factors.C, Sol_left.factors.D}; % tensor decompositition with temporal components, spectral components, spatial components, and subject components

% right hemisphere 
Sol_right=non_negative_cpd(tensor_right,r_components_r); % runnning a non-negative CPD
U_right={Sol_right.factors.A,Sol_right.factors.B, Sol_right.factors.C, Sol_right.factors.D}; % tensor decompositition with temporal components, spectral components, spatial components, and subject components

% left hemisphere
res = cpdres(tensor_left,U_left); % the residual between the tensor and its CPD approximation defined by U_left
relerr_left= frob(res)/frob(tensor_left); % the relative error between the tensor and its CPD approximation 

% right hemisphere
res = cpdres(tensor_right,U_right); % the residual between the tensor and its CPD approximation defined by U_right
relerr_right = frob(res)/frob(tensor_right); % the relative error between the tensor and its CPD approximation 

% save variables 
save(fullfile(saveDir,'CPD_all'),'tensor_left','tensor_right','condition_ind', ...
    'Sol_left','U_left','Sol_right','U_right','relerr_left','relerr_right', ...
    'r_components_l', 'r_components_r'...
    ); 

% plotting temporal component, spatial component, spectral component, and
% subject component
tick_time=(window_size-overlap)/fs; 
tick_freq=nfft/fs;
CPD_plots(U_left,r_components_l,true,condition_ind,conditions_title,saveDir,tick_time,tick_freq);
CPD_plots(U_right,r_components_r,false,condition_ind,conditions_title,saveDir,tick_time,tick_freq);

% plotting the temporal component, spatial component, spectral component, and the mean of each condition 
CPD_plots_mean(U_left,r_components_l,true,condition_ind,conditions_title,saveDir,tick_time,tick_freq);
CPD_plots_mean(U_right,r_components_r,false,condition_ind,conditions_title,saveDir,tick_time,tick_freq);

end 

% perfroms a non-negative CPD
function Sol=non_negative_cpd(Tens,r_components)
    model = struct;
    model.variables.a = rand(size(Tens,1),r_components);
    model.variables.b = rand(size(Tens,2),r_components);
    model.variables.c = rand(size(Tens,3),r_components);
    model.variables.d = rand(size(Tens,4),r_components);
    model.factors.A = {'a', @struct_nonneg};
    model.factors.B = {'b', @struct_nonneg};
    model.factors.C = {'c', @struct_nonneg};
    model.factors.D = {'d', @struct_nonneg};
    model.factorizations.tensor.data = Tens;
    model.factorizations.tensor.cpd = {'A', 'B', 'C','D'};
    Sol = sdf_nls(model);
end
