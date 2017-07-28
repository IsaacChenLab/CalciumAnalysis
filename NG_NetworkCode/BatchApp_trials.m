%runs BatchApp at multiple different connectivity thresholds
%and on different files
%requires modifications to the BatchApp and FC_crosscorr2.m
%need to be in same directory as BatchApp and params.mat

thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95];

base_path = '/Users/Nick/Documents/ChenLab/ChenLabData/FC Raw Trace testing/All cells/';
file_name = ' threshold - raw - all cells';

load('params.mat')

for i = 1:length(thresholds)
    path = strcat(base_path, num2str(thresholds(i)), file_name);
    
    %fprintf("path is %s\n", path);
    
    params.FC_thresh = thresholds(i);
    save('params.mat', '-append');
    ChenNetworkBatch(5, path);
    
end

fprintf("\tDone with trials\n");