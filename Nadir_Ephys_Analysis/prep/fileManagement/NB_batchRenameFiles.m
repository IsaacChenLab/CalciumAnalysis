% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% June 22, 2016

function [] = batchRenameFiles(numRecordings)

for rec = 0:numRecordings-1
    
    if rec < 10
        recordingDirectory = strcat('Recording_000', num2str(rec));
        cd(recordingDirectory);
        NB_ephysRename;
        cd ..;
    
    elseif rec < 100
        recordingDirectory = strcat('Recording_00', num2str(rec));
        cd(recordingDirectory);
        NB_ephysRename;
        cd ..;
        
    elseif rec < 1000
        recordingDirectory = strcat('Recording_0', num2str(rec));
        cd(recordingDirectory);
        NB_ephysRename;
        cd ..;
    
    else
        recordingDirectory = strcat('Recording', num2str(rec));
        cd(recordingDirectory);
        NB_ephysRename;
        cd ..;
    end
end 
        
    