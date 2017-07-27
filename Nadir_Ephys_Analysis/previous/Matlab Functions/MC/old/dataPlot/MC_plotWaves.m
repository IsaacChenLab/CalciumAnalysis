function MC_plotWaves(dataMat,Time,Mapping)

[nSamples,nChannels]=size(dataMat);

if nargin<2
    Time=1:nSamples;
end
if nargin<3,
    Mapping=[1:nChannels]';
end

maxs=max(dataMat);
mins=min(dataMat);
rngs=range(dataMat);

q=linspace(9,4,24);
fnt=round(q(nChannels));
figure;
bottom=0.06;
left=0.06; 
width=0.88;

for i=1:nChannels,
    %hc=axes('position',[left bottom width 0.88/nChannels-0.1/nChannels]);
    subplot(size(Mapping,1),size(Mapping,2),Mapping(i)); hold on;
    plot(Time,dataMat(:,i));
    set(gca,'FontSize',fnt);
    if size(Mapping,2)==1 & i~=nChannels
        set(gca,'XTick',[]);
    end
    set(gca,'YLim',[mins(i)-rngs(i)/7 maxs(i)+rngs(i)/7]);
    bottom=0.06+i*0.9/nChannels;
end

return;

