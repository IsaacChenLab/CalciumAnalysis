function MC_plotSharp2(ind)


if ind==1,
dcat('c369');
fl=[31 32]; fls=[31 32];
groups={[13:16 9:12 5:8 1:4],[21:24 17:20],[25:36]};
ch=15;
else
dcat('c370');
fl=[21 22]; fls=[21 22];
groups={[13:16 9:12 5:8 1:4],[21:24 17:20],[25:36],[37:44]};
ch=ind;
end

cdc;
window=[-1000 1000];
[eMat]=MC_loadEvent('Sharp',fls,ch);
[dataStream,totalMS,Sampling_Hz]=MC_openFile(fls(1));
if length(fl)==1
    xx=find(fls==fl);
    if xx==1
        eMat=eMat(1:totalMS);
    else
        eMat=eMat(totalMS+1 : end);
    end
end
for k=1:length(groups),
    figure; d=1;
    for i=groups{k}
        data=[];
        [SpikeMat,SpikeGrade]=MC_loadSpikes(fl,i);
        if ~isnan(SpikeMat)
            for ii=1:length(SpikeGrade)
                [counts_,mean_,std_]=MC_digOnDig(SpikeMat(:,ii),eMat,window);
                ca{d}=subplot(4,4,d); 
                bar([window(1):window(end)],counts_,'k');
                set(gca,'xlim',window);
                d=d+1; 
            end
        end
    end
    drawnow;
end

return;