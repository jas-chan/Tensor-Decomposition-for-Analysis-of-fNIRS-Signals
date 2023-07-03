%Forms the tensor and decomposes the tensor into components with Tucker
%decomposition (TD) 


%% Data directory and hyperparameters
clc; clear all; %clear comman window and workspace
data_dir_new='.\Text Files\';%access fNIRS datas
conditions={'Hand function_text','Hand motion_text','Robot function_text','Robot motion_text'}; % names of the folders for each condition
conditions_title={'Human hand, function', 'Human hand, non-function', 'Mechanical claw, function', ...
    'Mechanical claw, non-function'}; % actual names of the conditions
fs=50;  %sampling frequency
signal_length=1351; %number samples for trials 

saveDir=fullfile('Results','lmlra_decomp');   %directory of the tensor composition and decomposition 
if ~exist(saveDir)       %creating new directory  
    mkdir(saveDir)
end


%% Tensor Construction
tensor=zeros(signal_length,20,0); %establishing the tensor size
condition_ind=zeros(0); %keep track of the condition number
for conditionI=1:4 % going through all of the conditions 
    tensor_temp=read_condition_tensor(data_dir_new,conditions{conditionI}); %read the text files
    tensor_temp = fillmissing(tensor_temp,'constant',0); %dealing with the missing data (NaN)
    tensor(:,:,size(tensor,3)+1:size(tensor,3)+size(tensor_temp,3))=tensor_temp; %construct the tensor 
    condition_ind(end+1:end+size(tensor_temp,3))=conditionI; %keeping track of which subject goes with which condition 
end

%Splitting the tensor into 2, one for each hemisphere
tensor_left=tensor(:,1:10,:);     %creating left tensor from channels 1 through 10
tensor_right=tensor(:,11:end,:);  %creating right tensor from channels 11 to the end 


%% Tensor Decomposition 
%choosing the size of the core tensor
size_core_left = mlrankest(tensor_left); %determining the size of core (temporal by spectral by spatial by subject) for left hemisphere
size_core_right = mlrankest(tensor_right); %determining the size of core (temporal by spectral by spatial by subject) for right hemisphere
size_core_left(3)= size(tensor_left,3); %not decomposing the subject mode. This number should be the same number of subject files 
size_core_right(3)= size(tensor_right,3); %not decomposing the subject mode. This number should be the same number of subject files 

%decomposing (computes an intial soulution with the disired mulitilinear
%rank
[U_left,Sol_left] = lmlra(tensor_left,size_core_left);
[U_right,Sol_right] = lmlra(tensor_right,size_core_right);

% plotting temporal component, spatial component, spectral component, and
% subject component
TD_plots(U_left,size_core_left,true,condition_ind,conditions_title,saveDir,fs);
TD_plots(U_right,size_core_right,false,condition_ind,conditions_title,saveDir,fs);

% plotting the temporal component, spatial component, spectral component, and the mean of each condition 
TD_plots_mean(U_left,size_core_left,true,condition_ind,conditions_title,saveDir,fs);
TD_plots_mean(U_right,size_core_right,false,condition_ind,conditions_title,saveDir,fs);

save(fullfile(saveDir,'LMLRA_all')); % save variables 



