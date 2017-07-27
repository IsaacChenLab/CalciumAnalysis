function [SpikeMatAll,SpikeMats]=MC_loadSpikesAll1(file,ch)

global CAT;
cdc;
SpikeMat=nan;
file=MC_fileNumber(file);
if nargin==1
    nch=MC_NumElecInFile(file(1));
    ch=1:nch;
end
   
sp=[];
for cch=ch,
    [s]=MC_loadSpikes(file,cch,-2);
    if ~isnan(s)
        sp=[sp s];
    end
end

SpikeMats=sp;
SpikeMatAll=sum(sp,2);

return;