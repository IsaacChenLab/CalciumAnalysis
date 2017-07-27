function MC_scriptEvent(file,eventType,p1,Time)
% function MC_scriptSpike(file,ch,Time)
% spike sorting by file and then channels, meaning saves a mat file for each file with the
% channels inside it.
% call MC_spikes for start of spike sorting procedure, and saves the
% results in a file named Spike_[file #].
% thre results are 3 variables name SpikeTimes_[channel #],
% SpikeShapes_[channel #], and Mat_[channel #]
% can be called either by a specific filename, channel number, or can go
% over all channels in a file, or all on channels in all files.

global CAT;
global DATE;
global DATA_DIR;

cdc;
if nargin<4
    Time=nan;
end

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);    
nfl=sprintf('%d',file);

events_file=sprintf('Event_%s_%s_%s.mat',CAT,DATE,nfl);

        
[t_EventTimes,t_EvMat,t_Names]=MC_events(file,eventType,p1,Time);
% add to previous events
s=sprintf('EventTimes_%d=t_EventTimes;',p1);
eval(s);
s=sprintf('EventNames_%d=t_Names;',p1);
eval(s);
s=sprintf('EventMat_%d=t_EvMat;',p1);
eval(s);         

cd(DATA_DIR);
cd DigData;

if exist(events_file,'file')
    s=sprintf('save %s Ev* -append;',events_file);
    eval(s);    
else
    s=sprintf('save %s Ev*;',events_file);
    eval(s);    
end    



return;


