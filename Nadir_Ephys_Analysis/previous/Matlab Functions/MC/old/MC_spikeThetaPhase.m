function [SpikePhase,SpPhases]=MC_spikeThetaPhase(file,ch,Time)

SpikePhase=nan;

if nargin<3
    Time=nan;
end

[SpikeMat,SpikeGrade]=MC_loadSpikes(file,ch);

if isnan(SpikeMat)
    return;
end

[d,x]=MC_getTheta(file(1),ch,Time);
[ph1]=MC_phaseExtract(d,'extrema');
SpikeMat=SpikeMat(1:size(d,1),:);
SpikePhase=[];
for j=1:size(SpikeMat,2)
    SpPhases{j}=full(SpikeMat(:,j)).*ph1;
    SpPhases{j}=find(SpPhases{j});
    SpikePhase=[SpikePhase; SpPhases{j}];
end

return;



