%Reads channel indices from the .txt file
function data=load_data(file_path,channelNum) % input: the file path of where the .txt file is and the channels numbers to pull 
    %channelNum=[1:10] for left hemisphere and 11:20 for right.
    HbOInd=[10,4,1,13,22,7,16, 28, 25, 19, 40, 34, 31, 43, 52, 37, 46, 58, 55, 49]; %the column number of the channels of interest 
    T = readtable(file_path,'Delimiter',' ','ReadVariableNames',false,'HeaderLines',1, ...
    'Format',repmat('%f',[1,60])); % reads the.txt file
    data=table2array(T(:,HbOInd(channelNum)));
end