function [SpikeMatAll,SpikeMats,channels,grades,Inxs]=MC_loadSpikesAll(file,ch,spind)

global CAT;
if nargin<3
    spind=0;
end
SpikeMat=nan;
file=MC_fileNumber(file);
if nargin==1
    nch=MC_NumElecInFile(file(1));
    ch=1:nch;
end
   
sp=[]; channels=[]; grades=[]; Inxs=[];
for cch=ch,
    
    [s,g,ix]=MC_loadSpikes(file,cch,spind);
    if ~isempty(s)
        sp=[sp s];
        channels=[channels ones(1,size(s,2))*cch];
        grades=[grades g];
        Inxs=[Inxs ix];
    end
end

SpikeMats=sp;
SpikeMatAll=sum(sp,2);

return;