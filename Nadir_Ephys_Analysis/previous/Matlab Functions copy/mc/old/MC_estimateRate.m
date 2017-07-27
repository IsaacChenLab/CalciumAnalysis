function [Rate_Final,Time_Final]=MC_estimateRate(spikeMat,Bins)

if nargin<2
    Bins=50;
end

spikeMat=spikeMat(1:Bins*floor(size(spikeMat,1)/Bins),:);

nMs=size(spikeMat,1);
nSeconds=nMs/1000;
nCells=size(spikeMat,2);

spikeForSum=reshape(spikeMat,Bins,nMs*nCells/Bins);
spikeSum=sum(spikeForSum);
spikeMat=reshape(spikeSum,nMs/Bins,nCells);
Time=[1:size(spikeMat,1)]*Bins;

Rate=1e3*spikeMat/Bins;
Rate=full(Rate);
win=normpdf(-100:100,0,30); win=win/sum(win);
Rate=filtfilt(win,1,Rate);

for i=1:size(Rate,2),
    Rate_Final(:,i)=spline(Time,Rate(:,i),1:Time(end));
end
Time_Final=1:Time(end);

return;





