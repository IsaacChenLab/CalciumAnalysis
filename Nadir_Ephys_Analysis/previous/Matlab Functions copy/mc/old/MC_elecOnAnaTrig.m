function MC_elecOnAnaTrig(elecNum,anaCh,fileName,window,Time,h1,h2)
% function MC_elecOnAnaTrig(elecNum,anaCh,fileName,window,Time)

if nargin<4
    window=[-500 500];
end
if length(window)==1,
    window=[-round(window/2) round(window/2)];
end
if nargin<5
    Time=NaN;
end
[anaTimes anaInxs]=MC_getTriggerAnalogTimes(fileName,anaCh,Time,10,'prd');

for i=1:length(elecNum)
    [elecData,xtime,Hz]=MC_getElectrodes(fileName,elecNum(i),Time);
    window=window*Hz/1000;
    window=[window(1):window(end)];

    d=1; d1=[];
    for j=1:length(anaInxs),
        if (anaInxs(j) +window(1)) > 0 & (anaInxs(j)+window(end)) < size(elecData,1)
            d1(d,:)=elecData(anaInxs(j)+window,i)';
            d=d+1
        end
    end
    % normalize each trace to its baseline before trigger
    dm=mean(d1(:,1:abs(window(1)-2))')';
    d1=d1-repmat(dm,1,size(d1,2));
    dmean=nanmean(d1);
    dstd=nanstd(d1);
    dsem=dstd/sqrt(size(d1,1));
    
    %subplot(ceil(sqrt(length(elecNum))),ceil(sqrt(length(elecNum))),i); hold on;
    xax=window/(Hz/1000);
    if exist('h2','var')
        axes(h2);
    else
        figure; subplot(1,2,1); hold on;
    end
    plot(xax,d1');
    set(gca,'xlim',[xax(1) xax(end)]);
        
    if exist('h1','var')
        axes(h1);
    else
        subplot(1,2,2); hold on;
    end
    plot(xax,dmean,'b');
    [tt,xm]=max(dmean); 
    [tt,xn]=min(dmean);
    errorbar([xax(xm) xax(xn)],[dmean(xm) dmean(xn)],[dsem(xm) dsem(xn)],'r.');
    for_plot=dmean([1:round(length(dmean)/2)-100 round(length(dmean)/2)+100:end]);
    set(gca,'ylim',[mean(for_plot)-range(for_plot)*0.6 mean(for_plot)+range(for_plot)*0.6 ]);
    set(gca,'xlim',[xax(1) xax(end)]);
        
end
        
return;

    
