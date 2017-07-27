function MC_elecOnDigital(elecNum,digBitNum,fileName,window,Time)

if nargin<4
    window=[-500 500];
end
if length(window)==1,
    window=[-round(window/2) round(window/2)];
end
window=window*25;
window=[window(1):window(end)];
if nargin<5
    Time=NaN;
end

[digTimes digInxs]=MC_getDigitalTimes(fileName,digBitNum,Time);

for i=1:length(elecNum)
    [elecData,xtime]=MC_getElectrodes(fileName,elecNum(i),Time);
    d=1; d1=[];
    for j=1:length(digInxs),
        if (digInxs(j) +window(1)) > 0 & (digInxs(j)+window(end)) < size(elecData,1)
            d1(d,:)=elecData(digInxs(j)+window,i)';
            d=d+1;
        end
    end
    dmean=nanmean(d1);
    dstd=nanstd(d1);
    
    subplot(ceil(sqrt(length(elecNum))),ceil(sqrt(length(elecNum))),i); hold on;
    plot(window/25,dmean,'k');
    [tt,xm]=max(dmean); 
    [tt,xn]=min(dmean);
    errorbar([window(xm)/25 window(xn)/25],[dmean(xm) dmean(xn)],[dstd(xm) dstd(xn)],'k.');
        
end
        
return;

    
