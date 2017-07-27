function MC_spStart(file,ch,fileView,Time)
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

if nargin<4
    Time=nan;
end
if nargin<3
    fileView=file;
end
t=0;
[fileName,file]=MC_digFileName(file);
if ~exist(fileName,'file')
    eval(['save ' fileName ' t;']);
end


[spikeShapes_file,file]=MC_spShapeFileName(file);
if ~exist(spikeShapes_file,'file')
    eval(['save ' spikeShapes_file ' t;']);
end


for i=1:length(ch),
    s=sprintf('[SpikeTimes_%d,SpikeShapes_%d,SpikeMat_%d,SpikeGrades_%d,FileTimes]=MC_spSort(file,ch(i),fileView,Time);',ch(i),ch(i),ch(i),ch(i));
    eval(s);
    sss=eval(sprintf('SpikeMat_%d',ch(i)));
    if issparse(sss)
        s=sprintf('save %s -append SpikeT* SpikeG* SpikeM* FileTimes;',fileName);
        eval(s);    
        s=sprintf('save %s -append SpikeS*;',spikeShapes_file);
        eval(s);    
    end
    clear SpikeT* SpikeG* SpikeM* SpikeS*;
end

return;


