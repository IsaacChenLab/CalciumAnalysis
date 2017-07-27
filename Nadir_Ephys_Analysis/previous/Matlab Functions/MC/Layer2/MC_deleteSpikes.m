function [success]=MC_deleteSpikes(file,ch,sp)

% sp==-1 chooses the best spike (first by grade then by number of spikes)

global CAT;
global DATE;
global DATA_DIR;

success=0;

SpikeMat=[];
SpikeGrade=[];
SpikeTimes=[];

cd(DATA_DIR);
cd DigData;

fileName=MC_fileName(file);
file=MC_fileNumber(fileName);  
if length(file)==1 | ~isempty(find(diff(file)>1))
    nfl=sprintf('%d',file);
else
    nfl=sprintf('%d-%d',file(1),file(end));
end
spikes_file=sprintf('Spikes_%s_%s_%s.mat',CAT,DATE,nfl);
spikesShapes_file=sprintf('SpikesShapes_%s_%s_%s.mat',CAT,DATE,nfl);

s1=(sprintf('SpikeMat_%d',ch));
s2=(sprintf('SpikeGrades_%d',ch));
s3=(sprintf('SpikeTimes_%d',ch));
s4=(sprintf('SpikeShapes_%d',ch));

clear Sp* FileTimes;

cd(DATA_DIR);
cd DigData;
warning off all;
s=sprintf('load %s;',spikes_file);
eval(s);
s=sprintf('load %s;',spikesShapes_file);
eval(s);
warning on all;

if ~exist(s1,'var') | isempty(eval(s1)) | isnan(eval(s1))
    warning('Spike not found, doing nothing');
    return;
end

ns=size(eval(s1),2);
if nargin<3 | sp==-1 | ns==1
    s=sprintf('clear %s %s %s %s;',s1,s2,s3,s4);
    eval(s);
else
    inx=setdiff(1:nx,sp);
    s=sprintf('%s=%s(inx); %s=%s(inx); %s=%s(inx); %s=%s(inx);',s1,s1,s2,s2,s3,s3,s4,s4);
    eval(s);
end
s=sprintf('save %s SpikeT* SpikeG* SpikeM* FileTimes;',spikes_file);
eval(s);    
s=sprintf('save %s SpikeS*;',spikesShapes_file);
eval(s);  

success=1;

return;