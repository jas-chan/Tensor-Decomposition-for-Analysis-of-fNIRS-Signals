% Reads all the .txt files for one condition that is in one
% folder.
function tensor=read_condition_tensor(data_dir,condition_name)

    if exist(fullfile(strcat('conditions_data_',data_dir(end-2:end-1)),strcat(condition_name,'.mat')), 'file') == 2
         load (fullfile(strcat('conditions_data_',data_dir(end-2:end-1)),strcat(condition_name,'.mat')),'tensor');
         % If it already exists in .mat format, load the data. 
    else
        tensor=zeros(1351,20,1); % establishing the size: the number of time points by number of channels by 1 subject) 
        filesList= dir(fullfile(data_dir,condition_name));
        k=1;
        subj_files={};
        for i=1:size(filesList,1)
            if(filesList(i).isdir==0) %means file not folder
                subj_files(k)=cellstr(filesList(i).name);
                k=k+1;
            end
        end

        for fileI=1:length(subj_files)
            %load the data 
            data=load_data(fullfile(data_dir,condition_name,subj_files{fileI}),[1:20]); % input: path of .txt file and the number of channels pulled 
            tensor(:,:,fileI)=data*1000000; % convert to micromoles       
        end
        if ~exist(fullfile(strcat('conditions_data_',data_dir(end-2:end-1))))
            mkdir(fullfile(strcat('conditions_data_',data_dir(end-2:end-1))))
        end
        save(fullfile(strcat('conditions_data_',data_dir(end-2:end-1)),condition_name),'tensor','subj_files');
    end
end