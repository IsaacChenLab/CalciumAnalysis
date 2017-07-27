% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

function [] = NB_createFolders(numRecordings)

%
% [parentFolder deepestFolder] = fileparts(pwd);

% create folders for each separate recording
for i = 0:numRecordings-1;
    
    if i == 0;
        if ~exist('Recording_0000')
            mkdir('Recording_0000');
        end
    elseif i < 10;
        if ~exist(strcat('Recording_000', num2str(i)));
           mkdir(strcat('Recording_000', num2str(i)));
        end
    elseif i < 100;
        if ~exist(strcat('Recording_00', num2str(i)));
           mkdir(strcat('Recording_00', num2str(i)));
        end
    elseif i < 1000;
        if ~exist(strcat('Recording_0', num2str(i)));
            mkdir(strcat('Recording_0', num2str(i)));
        end
    else
        if ~exist(strcat('Recording', num2str(i)));
            mkdir(strcat('Recording', num2str(i)));
        end
    end
end